
function [Mat_state,stiff_mtx,Int_var]=...
    initial_constitutive(MAT_POINT,Mat_state,Int_var)

    global GEOMETRY SOLVER VARIABLE MATERIAL
    
    MODEL=MATERIAL.MODEL;
    Mat=MATERIAL.e;
    MAT=MATERIAL.MAT;
    
    df=GEOMETRY.df;
    dimf=GEOMETRY.f_dim;
    dims=GEOMETRY.s_dim;
    
    
    Kt=1;

    stiff_mtx = zeros(df*GEOMETRY.nodes);
    
    f_v       = zeros(dimf,1);
    f_old     = zeros(dimf,1);
    sig_0     = zeros(dims,1);
    
   
    TOL=1e-3;
    
    for e=1:GEOMETRY.mat_points
        
        tens_1=0;
        tens_2=0;
        
        if SOLVER.UW
            n=1-(1-MAT(16,Mat(e)))/MAT_POINT(e).J;
            dens=n*VARIABLE.rho_w+(1-n)*MAT(3,Mat(e));
        else
            dens=MAT(3,Mat(e))/MAT_POINT(e).J;
        end
        
        K0=MAT(26,Mat(e));
        if K0==0
            K0=1;
        end
        OCR=MAT(24,Mat(e));
        
        %% Initial Stress
        if SOLVER.INITIAL_COND(1)==1
            H=max(MAT_POINT(:).xg(2));
            tens_1=dens*VARIABLE.g*(H-MAT_POINT(e).xg(2));
            if SOLVER.UW==1
                Mat_state.pw(e,1)=VARIABLE.rho_w*VARIABLE.g*...
                    (H-MAT_POINT(e).xg(2));
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
        for i=1:dimf
            f_v(i,1)=Mat_state.F((e-1)*dimf + i,1);
            f_old(i,1)=Mat_state.F((e-1)*dimf + i,2);
        end           
        [F]=AUX.v2m(f_v);
        [Fold]=AUX.v2m(f_old);
        jacobians=det(F);
        
        if MODEL(Mat(e))>=2
            Sy=min(Int_var.Sy(e,2),OCR*press);
            Sy_r=Sy;
            Gamma=Int_var.gamma(e,2);
            dgamma=Int_var.dgamma(e,2);  
        end
             
        %% Calculate the deformation gradient state
        error(1)=1e32;
        ee=zeros(dims,1);
        e_fin=zeros(dims,1);
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
                [sig]=AUX.E2e(T);
                ds=sig-sig_0;
                
                D=A;

                if iter==1
                    ee(:,iter)=-A\ds;
                else
                    ee(:,iter)=ee(:,iter-1)-A\ds;
                end
                [E]=AUX.e2E(ee(:,iter));
            else
                D=A(1:3,1:3);

                [sig]=AUX.E2e(T);
                ds=sig(1:3)-sig_0(1:3);

                if iter==1
                    ee(1:3,iter)=-D\ds;
                else
                    ee(1:3,iter)=ee(1:3,iter-1)-D\ds;
                end
                [E]=AUX.e2E(ee(:,iter));
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
%         [E]=AUX.e2E(e_fin);
%         for i=1:3
%             F(i,i)=1/sqrt(exp(2*E(i,i)));
%         end
%         jacobians=det(F);
        
        
        %% Stiffness matrix
        [stiff_mtx]=stiff_mat(MAT_POINT,Mat_state,e,stiff_mtx,T,A);
        
        
        %% Store vectors
        T_vec(1,1)=T(1,1);
        T_vec(2,1)=T(2,2);
        T_vec(3,1)=T(3,3);
        T_vec(4,1)=T(1,2);
        for i=1:dims
            Mat_state.Sigma((e-1)*dims+i,1)=T_vec(i,1);
        end
        
        MAT_POINT(e).J    =   jacobians;
        if MODEL(Mat(e))>2
            Int_var.Sy(e,2)   =   Int_var.Sy(e,1);
            Int_var.Sy_r(e,2) =   Int_var.Sy_r(e,1);
            Int_var.P0(e,1)   =   press;
        end
        
        [be]=AUX.m2v(Be);
        [f]=AUX.m2v(F);
        for i=1:dimf
            Mat_state.Be((e-1)*dimf+i,2)=be(i,1);
            Mat_state.F((e-1)*dimf+i,2)=f(i,1);
        end 
        
%         if UW==1
%             [f_w]=AUX.m2v(F_w);
%             for i=1:dimf
%                 Mat_state.Fw((e-1)*dimf+i,2)=f_w(i,1);
%             end 
%         end
             
    end

end
