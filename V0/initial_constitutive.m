
function [Mat_state,stiff_mtx,Int_var,MAT_POINT]=...
    initial_constitutive(MAT_POINT,Mat_state,Int_var)

    global GEOMETRY SOLVER VARIABLE MATERIAL
    
    MODEL=MATERIAL.MODEL;
    Mat=MATERIAL.e;
    MAT=MATERIAL.MAT;
    
    df=GEOMETRY.df;
    dimf=GEOMETRY.f_dim;
    dims=GEOMETRY.s_dim;
        
    Kt=4;

    stiff_mtx = zeros(df*GEOMETRY.nodes);
    
    f_v       = zeros(dimf,1);
    f_old     = zeros(dimf,1);
    sig_0     = zeros(dims,1);
    
    
    for e=1:GEOMETRY.mat_points
        
        %% Initial Stress
        if SOLVER.INITIAL_COND(1)==1
            if SOLVER.UW
                n=1-(1-MAT(16,Mat(e)))/MAT_POINT(e).J;
                dens=n*VARIABLE.rho_w+(1-n)*MAT(3,Mat(e));
            else
                dens=MAT(3,Mat(e))/MAT_POINT(e).J;
            end
            
            H=max(MAT_POINT(:).xg(2));
            tens_1=dens*VARIABLE.g*(H-MAT_POINT(e).xg(2));
            if SOLVER.UW==1
                Mat_state.pw(e,1)=VARIABLE.rho_w*VARIABLE.g*...
                    (H-MAT_POINT(e).xg(2));
            end
        else
            tens_1=0;
        end
        tens_2=MAT(25,Mat(e));
        
%         if INITIAL_COND(3)~=0
%             Mat_state.pw(e,1)=INITIAL_COND(3);
%         end
        

        % MATRIX B and F
        b1=exp(MAT(23,Mat(e))/3)^2;
        Be=[b1 1e-15 0; 1e-15 b1 0; 0 0 b1];
        
        for i=1:dimf
            f_v(i,1)=Mat_state.F((e-1)*dimf + i,1);
            f_old(i,1)=Mat_state.F((e-1)*dimf + i,2);
        end           
        [F]=LIB.v2m(f_v);
        [Fold]=LIB.v2m(f_old);
        jacobians=det(F);
        
        
        % PRESSURE
        K0=MAT(26,Mat(e));
        if K0==0
            K0=1;
        end
        OCR=MAT(24,Mat(e));
        
        press=(1+2*K0)/3*tens_1+tens_2;
        
        
        if abs(press)>0
            
            SOLVER.step0=1;
            
            if MODEL(Mat(e))>=3
                nu=0.5;
                %nu=poisson(press,MODEL(Mat(e)),MAT(:,Mat(e)));
            else
                nu=MAT(2,Mat(e));
            end

            if SOLVER.AXI
                sig_0(2)=3*press/(1+2*K0);
                sig_0(1)=K0*sig_0(2);
                sig_0(4)=K0*sig_0(2);
            else
                sig_0(2)=3*press/(1+2*K0)*3/2/(1+nu);
                sig_0(1)=K0*sig_0(2);
                sig_0(4)=nu*K0*(sig_0(2)+sig_0(1));
            end


            if MODEL(Mat(e))>=2
                Sy=min(Int_var.Sy(e,2),OCR*press);
                Sy_r=Sy;
                Gamma=Int_var.gamma(e,2);
                dgamma=Int_var.dgamma(e,2); 
            end

            %% Calculate the deformation gradient state
            
            ee=zeros(dims,1);
            %e_fin=zeros(dims,1);
            iter=0;
            TOL=1e-3;
            imax=200;
            a=1;
            r0=abs(norm(sig_0));
            error=0;
            while iter <imax

                iter=iter+1;
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
                    elseif MODEL(Mat(e))>=4 && MODEL(Mat(e))<5
                        H=MAT(37,Mat(e));
                        [A,T,Gamma,dgamma,Sy,etaB,H,Be]=...
                            PZ(Kt,1,e,Gamma,dgamma,Sy,0,H,F,Fold,Be,press);
                    end
                    Int_var.Sy(e,1)     = Sy;
                    Int_var.gamma(e,1)  = Gamma;
                    Int_var.dgamma(e,1) = dgamma;
                end

                [sig]=LIB.E2e_in(T);
                ds=sig-sig_0;

                [CONVER,error,a,iter]=...
                    LIB.convergence(ds,r0,error,TOL,iter,imax,a);
                
                if SOLVER.AXI
                    if iter==1
                        ee(:,iter)=-a*A\ds;
                    else
                        ee(:,iter)=ee(:,iter-1)-a*A\ds;
                    end
                    [E]=LIB.e2E_in(ee(:,iter));
                else
                    D=A(1:3,1:3);
                    if iter==1
                        ee(1:3,iter)=-a*D\ds(1:3);
                    else
                        ee(1:3,iter)=ee(1:3,iter-1)-a*D\ds(1:3);
                    end
                    [E]=LIB.e2E_in(ee(:,iter));
                end

                for i=1:3
                    F(i,i)=exp(E(i,i));
                end
                jacobians=det(F);
                
                if CONVER==1     
                    break
                end

            end
           
            

            %% Stiffness matrix
            [stiff_mtx]=stiff_mat(MAT_POINT,Mat_state,e,stiff_mtx,T,A);

            %% Store vectors
            [T_vec]=LIB.E2e(T);
            for i=1:dims
                Mat_state.Sigma((e-1)*dims+i,1)=T_vec(i,1);
            end
            
            MAT_POINT(e).J    =   jacobians;
            if MODEL(Mat(e))>=2
                Int_var.Sy(e,2)   =   Int_var.Sy(e,1);
                Int_var.Sy_r(e,2) =   Int_var.Sy_r(e,1);
                Int_var.P0(e,1)   =   press;
                if MODEL(Mat(e))>=4
                    Int_var.H(e,1) =   H;
                    Int_var.H(e,2) =   H;
                    Int_var.eta(e,1) = etaB;
                    Int_var.eta(e,2) = etaB;
                end
            end
            
        
        else
            E=zeros(3);
            SOLVER.step0=0;
            
        end
            
        [be]=LIB.m2v(Be);
        [f]=LIB.m2v(F);
        for i=1:dimf
            Mat_state.Be((e-1)*dimf+i,2)=be(i,1);
            Mat_state.F((e-1)*dimf+i,2)=f(i,1);
        end 

            if SOLVER.UW
                Pw=SOLVER.INITIAL_COND(2);
                Mat_state.pw(:,1)=Pw;
                % Q calculation
                n=1-(1-MAT(16,Mat(e)))/jacobians;
                K_w=MAT(28,Mat(e));
                K_s=MAT(27,Mat(e));
                Q=1/(n/K_w+(1-n)/K_s);
                tr_e=E(1,1)+E(2,2)+E(3,3);
                tr_ew=-Pw/Q-tr_e;
                F_w=eye(3);
                if SOLVER.AXI
                    for i=1:3
                        F_w(i,i)=exp(tr_ew/3)^2;
                    end
                else
                    for i=1:2
                        F_w(i,i)=exp(tr_ew/2);
                    end
                end
                [f_w]=LIB.m2v(F_w);
                for i=1:dimf
                    Mat_state.Fw((e-1)*dimf+i,2)=f_w(i,1);
                end 
            end


    end

end

function nu=poisson(p,MODEL,MAT)

    if MODEL>=3 && MODEL<4 %CAMCLAY
        
        mu0   = MAT(4);   %mu0
        alfa  = MAT(20);  %alfa
        kappa = MAT(22);  %kappa
        
        K = -p/kappa;
        G = mu0-alfa*p;
        
        
    else %PZ

        khar  = MAT(29);
        ghar  = MAT(4);

        K=-khar*p;
        G=-ghar*p;
        
    end
    
    nu=(3*K-2*G)/2/(3*K+G);
    nu=min(0.5,max(nu,0));
end


