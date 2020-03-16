 classdef Constitutive
    methods(Static)
        function [stiff_mtx,Int_var,Mat_state]=...
            update(Kt,STEP,Int_var,Mat_state,MAT_POINT)

            global GEOMETRY SOLVER

            df=GEOMETRY.df;
            stiff_mtx = zeros(df*GEOMETRY.nodes);

            if SOLVER.FAIL==0
                for e=1:GEOMETRY.mat_points
                    % Constitutive calculation
                    [stiff_mtx,Int_var,Mat_state]=Constitutive.strain2const(...
                        e,STEP,Mat_state,Int_var,stiff_mtx,Kt,MAT_POINT);
                end
                
                % ----------------------------
                % Internal forces
                % ----------------------------
                [Mat_state]=Constitutive.internal_forces(MAT_POINT,Mat_state,STEP.BLCK);
            end   
        end

        function [stiff_mtx,Int_var,Mat_state]=strain2const(...
            e,STEP,Mat_state,Int_var,stiff_mtx,Kt,MAT_POINT)

            global GEOMETRY SOLVER MATERIAL

            ste=STEP.ste;
            BLCK=STEP.BLCK;

            MODEL=MATERIAL(BLCK).MODEL;
            Mat=GEOMETRY.material;
            MAT=MATERIAL(STEP.BLCK).MAT;

            dims=GEOMETRY.s_dim;

            % ----------------------------
            %% Strain preparation
            % ----------------------------
            if SOLVER.SMALL==0
                dimf=GEOMETRY.f_dim;

                f_v       = zeros(dimf,1);
                f_old     = zeros(dimf,1);
                be        = zeros(dimf,1); 
                %Vector F, Fp to Matrixes
                for i=1:dimf
                    f_v(i,1)=Mat_state.F((e-1)*dimf + i,1);
                    f_old(i,1)=Mat_state.F((e-1)*dimf + i,2);
                    be(i,1)=Mat_state.Be((e-1)*dimf + i,2);
                end           
                [F]=LIB.v2m(f_v);
                [Fold]=LIB.v2m(f_old);
                [Be_old]=LIB.v2m(be);
                J=det(F);

                if MODEL(Mat(e))<2
                    Be=F*F';
                end

                if MODEL(Mat(e))==0 || MODEL(Mat(e))>=2
                %%% Initial values  
                    if MODEL(Mat(e))>=2
                        % Predictor tensors
                        Fincr=F/Fold;
                        % Compute Trial left cauchy-Green
                        Be = Fincr*Be_old*Fincr';   
                        if isnan(BeTr)
                            error('Error in Green-Lagrange tensor of elem e %i \n',e);
                        end
                        if MODEL(Mat(e))>=4
                            E_tot = logm(F*F')/2;
                            E_0   = logm(Fold*Fold')/2;
                            deps =  E_tot- E_0;
                        end
                    end
                    % Elastic trial strain
                    Ee = logm(Be)/2;
                    if isnan(Ee)
                        error('Error in small strain tensor of elem e %i \n',e);
                    elseif isreal(Ee)==0
                        error('Complex in small strain tensor of elem e %i \n',e);
                    end
                end
            else
                e_v       = zeros(dims,1);
                e_old     = zeros(dims,1);
                e_e       = zeros(dims,1); 
                %Vector F, Fp to Matrixes
                for i=1:dims
                    e_v(i,1)  = Mat_state.Es((e-1)*dims + i,1);
                    e_old(i,1)= Mat_state.Es((e-1)*dims + i,2);
                    e_e(i,1)  = Mat_state.Es_e((e-1)*dims + i,2);
                end 
                [Et]=LIB.e2E(e_v);
                [Et_old]=LIB.e2E(e_old);
                [Ee_old]=LIB.e2E(e_e);

                deps=Et-Et_old;
                Ee=Ee_old+deps;
            end

            % ----------------------------
            %% Constitutive Calculation
            % ----------------------------
            if MODEL(Mat(e))<2
                if MODEL(Mat(e))==0
                    [A,T]=Saint_Venant(Kt,e,Ee,BLCK);
                elseif MODEL(Mat(e))<2 && MODEL(Mat(e))>=1
                    if SOLVER.SMALL==0
                        [A,T]=Neo_Hookean(Kt,e,Be,J,BLCK);
                    else
                        [A,T]=Saint_Venant(Kt,e,Ee,BLCK);
                    end
                end
            else        
                Sy      = Int_var.Sy(e,2);
                Gamma   = Int_var.gamma(e,2);
                dgamma  = Int_var.dgamma(e,2);

                if MODEL(Mat(e))>=2 && MODEL(Mat(e))<3
                    [A,T,Gamma,dgamma,Sy,Ee]=...
                        Drucker_prager(Kt,e,Gamma,dgamma,Sy,Ee,BLCK);
                elseif MODEL(Mat(e))>=3 && MODEL(Mat(e))<4
                    %P0      = Int_var.P0(e);
                    Sy_r    = Int_var.Sy_r(e,2);
                    [A,T,Gamma,dgamma,Sy,Sy_r,Ee]=...
                        M_Cam_Clay(Kt,ste,e,Gamma,dgamma,Sy,Sy_r,Ee,BLCK);
                    Int_var.Sy_r(e,1) = Sy_r;
                elseif MODEL(Mat(e))>=4 && MODEL(Mat(e))<5
                    %P0      = Int_var.P0(e);
                    H       = Int_var.H(e,2);
                    etaB    = Int_var.eta(e,2);
                    epsvol  = Int_var.epsv(e,2);
                    [A,T,Gamma,epsvol,dgamma,Sy,etaB,H,Ee]=...
                            PZ(Kt,ste,e,Gamma,epsvol,dgamma,Sy,etaB,H,Ee,deps,BLCK);
                    Int_var.epsv(e,1)= epsvol;
                    Int_var.H(e,1)   = H;
                    Int_var.eta(e,1) = etaB;
                end
                Int_var.Sy(e,1)     = Sy;
                Int_var.gamma(e,1)  = Gamma;
                Int_var.dgamma(e,1) = dgamma;
            end

            if SOLVER.SMALL==0 && MODEL(Mat(e))>=2
                Be = expm(2*Ee); 
                T = T/J;        
            end

            AA=isreal(T);
            if  AA(1)==0
                fprintf('Error with stress, in Constitutive file, element %i \n',e);
                SOLVER.FAIL=1;
            else
            % ----------------------------
            %% Stress storage
            % ----------------------------     
                [sig]=LIB.E2e(T);
                for i=1:dims
                    Mat_state.Sigma((e-1)*dims+i,1)=sig(i,1);
                end
                if SOLVER.SMALL==0
                    [be]=LIB.m2v(Be);
                    for i=1:dimf
                        Mat_state.Be((e-1)*dimf+i,1)=be(i,1);
                    end 
                else
                    [evec]=LIB.E2e(Ee);
                    for i=1:dims
                        Mat_state.Es_e((e-1)*dims+i,1)=evec(i,1);
                    end
                end

            % ----------------------------
            %% Stiffness matrix
            % ----------------------------
                if Kt==1 || Kt==2 || Kt==4
                    [stiff_mtx]=...
                        stiff_mat(MAT_POINT,Mat_state,e,stiff_mtx,T,A,BLCK);
                end           
            end

            if SOLVER.UW==1
            % ----------------------------
            %% Pore Water Pressure Calculation matrix
            % ----------------------------
                % Q calculation
                n=1-(1-MAT(16,Mat(e)))/MAT_POINT{1}(e).J;
                K_w=MAT(28,Mat(e));
                K_s=MAT(27,Mat(e));
                Q=1/(n/K_w+(1-n)/K_s);

                %%% Strains %%%
                if SOLVER.SMALL==0
                    f_v_w     = zeros(dimf,1);
                    for i=1:dimf
                        f_v_w(i,1)=Mat_state.Fw((e-1)*dimf + i,1);
                    end           
                    [F_w]=LIB.v2m(f_v_w);

                    bw=F_w*F_w';
                    b=F*F';
                    % Compute principal deformation and direction
                    [~,eigva_b] = eig(b);
                    tr_e = log(sqrt(eigva_b(1,1)))+log(sqrt(eigva_b(2,2)))+...
                            log(sqrt(eigva_b(3,3)));
                    [~,eigva_bw] = eig(bw);
                    tr_ew = log(sqrt(eigva_bw(1,1)))+log(sqrt(eigva_bw(2,2)))+...
                            log(sqrt(eigva_bw(3,3)));
                else
                    e_w       = zeros(dims,1); 
                    %Vector Ew to matrix
                    for i=1:dims
                        e_w(i,1)  = Mat_state.Esw((e-1)*dims + i,1);
                    end 
                    [Ew]=LIB.e2E(e_w);
                    tr_e  = trace(Et);
                    tr_ew = trace(Ew);
                end

                % Pore pressure
                Mat_state.pw(e,1)=-Q*(tr_e+tr_ew);
            end
        end

        function [Mat_state]=internal_forces(MAT_POINT,Mat_state,BLCK)

            global GEOMETRY SOLVER

            sig=zeros(4,1);
            Mat_state.fint(:,1) = zeros(GEOMETRY.nodes*GEOMETRY.df,1);
            sp=GEOMETRY.sp;
            df=GEOMETRY.df;

            for e=1:GEOMETRY.mat_points

                nd=MAT_POINT{1}(e).near;
                B_=MAT_POINT{1}(e).B;
                nn =length(nd);

                if SOLVER.UW>0
                    ndw=MAT_POINT{2}(e).near;
                    B_w=MAT_POINT{2}(e).B;
                    nnw =length(ndw); 
                end

                % Derivatives
                if SOLVER.AXI
                    sh=zeros(3,nn);
                    for i=1:nn
                        sh(1,i)=B_(1,i*2-1);
                        sh(2,i)=B_(2,i*2);
                        sh(3,i)=B_(4,i*2-1);
                    end
                    if SOLVER.UW>0
                        shw=zeros(3,nnw);
                        for i=1:nnw
                            shw(1,i)=B_w(1,i*2-1);
                            shw(2,i)=B_w(2,i*2);
                            shw(3,i)=B_w(4,i*2-1);
                        end
                    end
                else
                    sh=zeros(2,nn);
                    for i=1:nn
                        sh(1,i)=B_(1,i*2-1);
                        sh(2,i)=B_(2,i*2);
                    end
                    if SOLVER.UW>0
                        shw=zeros(2,nnw);
                        for i=1:nnw
                            shw(1,i)=B_w(1,i*2-1);
                            shw(2,i)=B_w(2,i*2);
                        end
                    end
                end

                % Stress
                for i=1:4
                    sig(i,1)=Mat_state.Sigma((e-1)*4+i,1)-Mat_state.Sigma((e-1)*4+i,3);
                end

                [T]=LIB.e2E(sig);

                % ----------------------------
                % Internal forces
                % ----------------------------
                volume=GEOMETRY.Area(e)*MAT_POINT{1}(e).J;
                if SOLVER.AXI
                    mat=[1 0 1; 0 1 0];
                    Tt=mat*T;
                    vol=2*pi*MAT_POINT{1}(e).xg(1)*volume;
                else
                    Tt=T(1:2,1:2);
                    vol=volume;
                end
                int_forces_1=Tt*sh*vol;

                if SOLVER.UW==1
                    sh2=sh(1:2,:);
                    sh2w=shw(1:2,:);
                    if SOLVER.AXI
                        sh2(1,:)=sh2(1,:)+sh(3,:);
                        sh2w(1,:)=sh2w(1,:)+shw(3,:);
                    end
                    int_forces_2=sh2*(Mat_state.pw(e,1)-Mat_state.pw(e,3))*vol;
                    int_forces_3=sh2w*(Mat_state.pw(e,1)-Mat_state.pw(e,3))*vol;
                    if SOLVER.IMPLICIT(BLCK)==0
                        for i=1:nn
                           nod=nd(i);
                           for j=1:sp
                                Mat_state.fint(nod*df+1-sp-j,1)=...
                                    Mat_state.fint(nod*df+1-sp-j,1)-int_forces_1(3-j,i);
                                Mat_state.fint(nod*df+1-j,1)=...
                                    Mat_state.fint(nod*df+1-j,1)-int_forces_2(3-j,i);
                           end
                        end
                    else
                        for i=1:nn
                           nod=nd(i);
                           for j=1:sp
                                Mat_state.fint(nod*df-1-j,1)=...
                                    Mat_state.fint(nod*df-1-j,1)-int_forces_1(3-j,i)+...
                                    int_forces_2(3-j,i);
                           end
                        end
                        for i=1:nnw
                           nod=ndw(i);
                           for j=1:sp
                                Mat_state.fint(nod*df+1-j,1)=...
                                    Mat_state.fint(nod*df+1-j,1)+int_forces_3(3-j,i);
                           end
                        end
                    end            
                elseif SOLVER.UW==2
                    if SOLVER.AXI
                        m=[1 1 0 1];
                    else
                        m=[1 1 0];
                    end
                    div=B_'*m';
                    dN=zeros(2,nnw);
                    for j=1:nnw
                        dN(1,j)=B_w(1,(j-1)*sp+1);
                        dN(2,j)=B_w(2,(j-1)*sp+2);
                    end
                    int_forces_2=div*(Mat_state.pw(e,1)-Mat_state.pw(e,3))*vol;
                    dPw=Mat_state.dpw((e-1)*sp+1:e*sp,1);
                    int_forces_3=Mat_state.k(e)*dN'*dPw*vol;
                    for i=1:nn
                       nod=nd(i);
                       for j=1:sp
                            Mat_state.fint(nod*df-j,1)=...
                                Mat_state.fint(nod*df-j,1)-int_forces_1(3-j,i)+...
                                int_forces_2(i*sp+1-j,1);
                       end
                    end
                    for i=1:nnw
                       nod=ndw(i);
                       Mat_state.fint(nod*df,1)=...
                                Mat_state.fint(nod*df,1)-int_forces_3(i,1);
                    end
                else
                    for i=1:nn
                       nod=nd(i);
                       for j=1:sp
                            Mat_state.fint(nod*sp+1-j,1)=Mat_state.fint(nod*sp+1-j,1)+...
                                int_forces_1(3-j,i);
                       end
                    end
                end
                clear sh
            end
        end
        
    end
 end