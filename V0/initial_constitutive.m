
function [Mat_state,stiff_mtx,Int_var]=...
    initial_constitutive(Shape_function,Mat_state,Int_var)

    global GEOMETRY SOLVER VARIABLE MATERIAL
    
    MODEL=MATERIAL.MODEL;
    Mat=MATERIAL.e;
    MAT=MATERIAL.MAT;
    
    
    Kt=1;

    stiff_mtx = zeros(GEOMETRY.df*GEOMETRY.nodes,GEOMETRY.df*GEOMETRY.nodes);
    
    f_v       = zeros(GEOMETRY.sp*GEOMETRY.sp+1,1);
    f_old     = zeros(GEOMETRY.sp*GEOMETRY.sp+1,1);
    sig_0     = zeros(4,1);
    
   
    TOL=1e-3;
    
    for e=1:GEOMETRY.elements
        
        tens_1=0;
        tens_2=0;
        
        if SOLVER.UW
            n=1-(1-MAT(16,Mat(e)))/Mat_state.J(e);
            dens=n*VARIABLE.rho_w+(1-n)*MAT(3,Mat(e));
        else
            dens=MAT(3,Mat(e))/Mat_state.J(e);
        end
        
        K0=MAT(26,Mat(e));
        if K0==0
            K0=1;
        end
        OCR=MAT(24,Mat(e));
        
        %% Initial Stress
        if SOLVER.INITIAL_COND(1)==1
            H=max(Mat_state.xg(:,2));
            tens_1=dens*VARIABLE.g*(H-Mat_state.xg(e,2));
            if SOLVER.UW==1
                Mat_state.pw(e,1)=VARIABLE.rho_w*VARIABLE.g*...
                    (H-Mat_state.xg(e,2));
            end
        end
        

        tens_2=MAT(25,Mat(e));
        
%         if INITIAL_COND(3)~=0
%             Mat_state.pw(e,1)=INITIAL_COND(3);
%         end
        
        b1=exp(MAT(23,Mat(e))/3)^2;
        Be=[b1 1e-15 0; 1e-15 b1 0; 0 0 b1];
        
        
        press=(1+2*K0)/3*tens_1+tens_2;
        
        sig_0(2)=3*press/(1+2*K0);
        sig_0(1)=K0*sig_0(2);
        sig_0(4)=K0*sig_0(2);
        
        
        %% Vector F, Fp to Matrixes
        for i=1:5
            f_v(i,1)=Mat_state.F((e-1)*5 + i,1);
            f_old(i,1)=Mat_state.F((e-1)*5 + i,2);
        end           
        [F]=v2m(f_v,GEOMETRY.sp);
        [Fold]=v2m(f_old,GEOMETRY.sp);
        jacobians=det(F);
        
        if MODEL(Mat(e))>2
            Sy=min(Int_var.Sy(e,2),OCR*press);
            Sy_r=Sy;
            Gamma=Int_var.gamma(e,2);
            dgamma=Int_var.dgamma(e,2);  
        end
             
        %% Calculate the deformation gradient state
        error(1)=1e32;
        ee=zeros(4,1);
        e_fin=zeros(4,1);
        iter=1;
        while error(iter) > TOL
        
            if MODEL(Mat(e))<2
                if MODEL(Mat(e))==0
                    [A,T,Be]=Saint_Venant(Kt,e,F);
                elseif MODEL(Mat(e))<2 && MODEL(Mat(e))>=1
                    [A,T,Be]=Neo_Hookean(Kt,e,F,jacobians);
                end
            else        
                if MODEL(Mat(e))>=2 && MODEL(Mat(e))<3
                    [A,T,Gamma,dgamma,Sy,Be]=...
                        Drucker_prager(Kt,e,Gamma,dgamma,Sy,F,Be,Fold);
                elseif MODEL(Mat(e))>=3 && MODEL(Mat(e))<4
                    [A,T,Gamma,dgamma,Sy,Sy_r,Be]=...
                        M_Cam_Clay(Kt,1,e,Gamma,dgamma,Sy,Sy_r,F,Fold,Be,press);
                    Int_var.Sy_r(e,1) = Sy_r;
                end
                Int_var.Sy(e,1)     = Sy;
                Int_var.gamma(e,1)  = Gamma;
                Int_var.dgamma(e,1) = dgamma;
            end
            
            if SOLVER.AXI
                [sig]=E2e(T);
                ds=sig-sig_0;
                
                D=A;

                if iter==1
                    ee(:,iter)=-A\ds;
                else
                    ee(:,iter)=ee(:,iter-1)-A\ds;
                end
                [E]=e2E(ee(:,iter));
            else
                D=A(1:3,1:3);

                [sig]=E2e(T);
                ds=sig(1:3)-sig_0(1:3);

                if iter==1
                    ee(1:3,iter)=-D\ds;
                else
                    ee(1:3,iter)=ee(1:3,iter-1)-D\ds;
                end
                [E]=e2E(ee(:,iter));
            end

            for i=1:3
                F(i,i)=sqrt(exp(2*E(i,i)));
            end
            jacobians=det(F);
            
            iter=iter+1;
            error(iter)=abs(norm(ds));

        end
        
        %% DEFORMATION GRADIENT
%         e_fin=D\sig;  
%         [E]=e2E(e_fin);
%         for i=1:3
%             F(i,i)=1/sqrt(exp(2*E(i,i)));
%         end
%         jacobians=det(F);
        
        
        %% Stiffness matrix
        [stiff_mtx]=stiff_mat(Shape_function,Mat_state,e,stiff_mtx,T,A);
        
        
        %% Store vectors
        T_vec(1,1)=T(1,1);
        T_vec(2,1)=T(2,2);
        T_vec(3,1)=T(3,3);
        T_vec(4,1)=T(1,2);
        for i=1:4
            Mat_state.Sigma((e-1)*4+i,1)=T_vec(i,1);
        end
        
        Mat_state.J(e)    =   jacobians;
        if MODEL(Mat(e))>2
            Int_var.Sy(e,2)   =   Int_var.Sy(e,1);
            Int_var.Sy_r(e,2) =   Int_var.Sy_r(e,1);
            Int_var.P0(e,1)   =   press;
        end
        
        [be]=m2v(Be,GEOMETRY.sp);
        [f]=m2v(F,GEOMETRY.sp);
        for i=1:5
            Mat_state.Be((e-1)*5+i,2)=be(i,1);
            Mat_state.F((e-1)*5+i,2)=f(i,1);
        end 
        
%         if UW==1
%             [f_w]=m2v(F_w,sp);
%             for i=1:5
%                 Mat_state.Fw((e-1)*5+i,2)=f_w(i,1);
%             end 
%         end
             
    end

end

function [E]=e2E(e)
   
    E=zeros(3,3);

    %Build matrix
    E(1,1)=e(1);
    E(2,2)=e(2);
    E(3,3)=e(4);
    E(1,2)=e(3);
    E(2,1)=e(3);
end

function [e]=E2e(E)
   
    e=zeros(4,1);

    %Build vector
    e(1)=E(1,1);
    e(2)=E(2,2);
    e(3)=E(1,2);
    e(4)=E(3,3);%+E(2,1);
end
