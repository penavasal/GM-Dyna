function [stiff_mtx,Int_var,Mat_state,FAIL]=...
    Constitutive(Kt,ste,Int_var,Mat_state,Shape_function,FAIL)
 
    global GEOMETRY SOLVER MATERIAL
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    
    Mat=MATERIAL.e;
    MODEL=MATERIAL.MODEL;
    MAT=MATERIAL.MAT;
    
    f_v       = zeros(sp*sp+1,1);
    f_old     = zeros(sp*sp+1,1);
    if SOLVER.UW
        f_v_w     = zeros(sp*sp+1,1);
    end
    be        = zeros(sp*sp+1,1); 
    stiff_mtx = zeros(df*GEOMETRY.nodes,df*GEOMETRY.nodes);

    if FAIL==0
        for e=1:GEOMETRY.elements
                        
            if SOLVER.UW
                for i=1:5
                    f_v_w(i,1)=Mat_state.Fw((e-1)*5 + i,1);
                end           
                [F_w]=v2m(f_v_w,sp);
            end
                
            %Vector F, Fp to Matrixes
            for i=1:5
                f_v(i,1)=Mat_state.F((e-1)*5 + i,1);
                f_old(i,1)=Mat_state.F((e-1)*5 + i,2);
                be(i,1)=Mat_state.Be((e-1)*5 + i,2);
            end           
            [F]=v2m(f_v,sp);
            [Fold]=v2m(f_old,sp);
            [Be]=v2m(be,sp);
            
            if MODEL(Mat(e))<2
                if MODEL(Mat(e))==0
                    [A,T,Be]=Saint_Venant(Kt,e,F);
                elseif MODEL(Mat(e))<2 && MODEL(Mat(e))>=1
                    [A,T,Be]=Neo_Hookean(Kt,e,F,Mat_state.J(e));
                end
            else        
                Sy      = Int_var.Sy(e,2);
                Gamma   = Int_var.gamma(e,2);
                dgamma  = Int_var.dgamma(e,2);
		
                if MODEL(Mat(e))>=2 && MODEL(Mat(e))<3
                    [A,T,Gamma,dgamma,Sy,Be]=...
                        Drucker_prager(Kt,e,Gamma,dgamma,Sy,F,Be,Fold);
                elseif MODEL(Mat(e))>=3 && MODEL(Mat(e))<4
                    P0      = Int_var.P0(e);
                    Sy_r     = Int_var.Sy_r(e,2);
                    [A,T,Gamma,dgamma,Sy,Sy_r,Be]=...
                        M_Cam_Clay(Kt,ste,e,Gamma,dgamma,Sy,Sy_r,F,Fold,Be,P0);
                    Int_var.Sy_r(e,1) = Sy_r;
                end
                Int_var.Sy(e,1)     = Sy;
                Int_var.gamma(e,1)  = Gamma;
                Int_var.dgamma(e,1) = dgamma;
            end

            
            if SOLVER.UW
                % Q calculation
                n=1-(1-MAT(16,Mat(e)))/Mat_state.J(e);
                K_w=MAT(28,Mat(e));
                K_s=MAT(27,Mat(e));
                Q=1/(n/K_w+(1-n)/K_s);
                % Pressure
                bw=F_w*F_w';
                b=F*F';
                % Compute principal deformation and direction
                [~,eigva_b] = eig(b);
                tr_e = log(sqrt(eigva_b(1,1)))+log(sqrt(eigva_b(2,2)))+...
                        log(sqrt(eigva_b(3,3)));
                [~,eigva_bw] = eig(bw);
                tr_ew = log(sqrt(eigva_bw(1,1)))+log(sqrt(eigva_bw(2,2)))+...
                        log(sqrt(eigva_bw(3,3)));
                % Pore pressure
                Mat_state.pw(e,1)=-Q*(tr_e+tr_ew);
            end

            AA=isreal(T);
            if  AA(1)==0
                FAIL=1;
                break;
            else
                [sig]=E2e(T);
                [be]=m2v(Be,sp);
                for i=1:5
                    Mat_state.Be((e-1)*5+i,1)=be(i,1);
                end 
                for i=1:4
                    Mat_state.Sigma((e-1)*4+i,1)=sig(i,1);
                end
                
                % ----------------------------
                % Stiffness matrix
                % ----------------------------
                if Kt==1 || Kt==2 || Kt==4
                    [stiff_mtx]=...
                        stiff_mat(Shape_function,Mat_state,e,stiff_mtx,T,A);
                end           
            end
        end
        % ----------------------------
        % Internal forces
        % ----------------------------
        [Mat_state]=internal_forces(Shape_function,Mat_state);
    end   
end

function [e]=E2e(E)
   
    e=zeros(4,1);

    %Build vector
    e(1)=E(1,1);
    e(2)=E(2,2);
    e(3)=E(3,3);
    e(4)=E(1,2);%+E(2,1);
end
