 classdef VECTORS
    methods(Static)

        % Store
        function GLOBAL=Store(GLOBAL,STEP,Mat_state,Disp_field,...
            Int_var,MAT_POINT,out_list1)

            global SOLVER GEOMETRY
            
            ste_p=STEP.ste_p;

            [out_list2]=VECTORS.reaction(Mat_state.fint,STEP);
            GLOBAL.OutputList(ste_p,:)=out_list1+out_list2;

            GLOBAL.d(:,ste_p)   = Disp_field.d(:,1);
            GLOBAL.a(:,ste_p)   = Disp_field.a(:,1);
            GLOBAL.v(:,ste_p)   = Disp_field.v(:,1);

            GLOBAL.gamma(:,ste_p)   = Int_var.gamma(:,1);
            GLOBAL.epsv(:,ste_p)    = Int_var.epsv(:,1);
            GLOBAL.dgamma(:,ste_p)  = Int_var.dgamma(:,1);
            GLOBAL.H(:,ste_p)       = Int_var.H(:,1);
            GLOBAL.eta(:,ste_p)     = Int_var.eta(:,1);
            GLOBAL.Sy(:,ste_p)      = Int_var.Sy(:,1);
            GLOBAL.Sy_r(:,ste_p)    = Int_var.Sy_r(:,1);
            GLOBAL.H(:,ste_p)       = Int_var.H(:,1);
            GLOBAL.eta(:,ste_p)     = Int_var.eta(:,1);

            GLOBAL.xg(:,ste_p)      = LIB.reshape_S2list(MAT_POINT{1},'xg');

            [GLOBAL.gamma_nds(:,ste_p)]=VECTORS.Ep2Ep_n...
                (GLOBAL.gamma,MAT_POINT,ste_p);

            for e=1:GEOMETRY.mat_points
                [GLOBAL.Ps(e,ste_p),GLOBAL.Qs(e,ste_p)]=...
                    VECTORS.invar(Mat_state.Sigma(:,1),e);   %PRESSURE
            end
            
            [GLOBAL.Es(:,ste_p),GLOBAL.Es_p(:,ste_p)]=...
                VECTORS.strains(Mat_state);
            [GLOBAL.J(:,ste_p)]=LIB.S2list(MAT_POINT{1},'J'); %JACOBIAN
            GLOBAL.Sigma(:,ste_p)   = Mat_state.Sigma(:,1);    %STRESS
            if SOLVER.SMALL==0 
                GLOBAL.F(:,ste_p)       = Mat_state.F(:,1);  %DEFORMATION GRADIENT
                GLOBAL.Be(:,ste_p)      = Mat_state.Be(:,1); %LEFT CAUCHY GREEN 
            end
            
            if SOLVER.UW==1
                GLOBAL.pw(:,ste_p) = Mat_state.pw(:,1);
                if SOLVER.SMALL==0
                    GLOBAL.Fw(:,ste_p) = Mat_state.Fw(:,1);
                else
                    GLOBAL.Esw(:,ste_p) = Mat_state.Esw(:,1);
                end
            elseif SOLVER.UW==2
                GLOBAL.pw(:,ste_p) = Mat_state.pw(:,1);
                GLOBAL.dpw(:,ste_p) = Mat_state.dpw(:,1);
            end
            
            if SOLVER.FRAC
                GLOBAL.w(:,ste_p)      = Mat_state.w(:,1);
                GLOBAL.status(:,ste_p) = Mat_state.status(:,1);
                if SOLVER.FRAC>1
                    GLOBAL.e_ini(:,ste_p) = Mat_state.e_ini(:,1);
                end
                GLOBAL.E.K(ste_p,SOLVER.BODIES)=0;
                GLOBAL.E.W(ste_p,SOLVER.BODIES)=0;
                GLOBAL.E.D(ste_p,1)=STEP.ENERGY.D;
                [GLOBAL.E]=FRAC.compute_energy(...
                    MAT_POINT,STEP,Mat_state,GLOBAL.E,Disp_field.v);
                
                GLOBAL.Force(ste_p,SOLVER.BODIES*GEOMETRY.sp)=0;
                [GLOBAL.Force]=FRAC.forces(MAT_POINT,STEP,GLOBAL.Force,Disp_field.v);
            end

            GLOBAL.tp(ste_p,1)=STEP.t;
            GLOBAL.ste_p=ste_p;

        end

        % Update
        function [Disp_field,Mat_state,Int_var]=Update(...
                Disp_field,Mat_state,Int_var)
            
            global SOLVER
            
            Disp_field.d(:,3)=Disp_field.d(:,2);
            Disp_field.d(:,2)=Disp_field.d(:,1);
            Disp_field.v(:,2)=Disp_field.v(:,1);
            Disp_field.a(:,2)=Disp_field.a(:,1);

            Mat_state.Sigma(:,2)=Mat_state.Sigma(:,1);                
            Mat_state.fint(:,2)=Mat_state.fint(:,1);
            
            if SOLVER.SMALL==0
                Mat_state.F(:,2)=Mat_state.F(:,1);
                Mat_state.Be(:,2)=Mat_state.Be(:,1);
            else
                Mat_state.Es(:,2)=Mat_state.Es(:,1);
                Mat_state.Es_e(:,2)=Mat_state.Es_e(:,1);
            end

            if SOLVER.UW==1
                Mat_state.pw(:,2)=Mat_state.pw(:,1);
                if SOLVER.SMALL==0
                    Mat_state.Fw(:,2)=Mat_state.Fw(:,1);
                else
                    Mat_state.Esw(:,2)=Mat_state.Esw(:,1);
                end
            elseif SOLVER.UW==2
                Mat_state.pw(:,2)=Mat_state.pw(:,1);
                Mat_state.dpw(:,2)=Mat_state.dpw(:,1);
            end
            
            if SOLVER.FRAC
                Mat_state.status(:,2)=Mat_state.status(:,1);
            end

            Int_var.gamma(:,2)  = Int_var.gamma(:,1);
            Int_var.epsv(:,2)   = Int_var.epsv(:,1);
            Int_var.H(:,2)      = Int_var.H(:,1);
            Int_var.eta(:,2)    = Int_var.eta(:,1);
            Int_var.Sy(:,2)     = Int_var.Sy(:,1);
            Int_var.Sy_r(:,2)   = Int_var.Sy_r(:,1);
            Int_var.dgamma(:,2) = Int_var.dgamma(:,1); 
        end
        
        % Update initial
        function [Disp_field,Mat_state,Int_var,stiff_mtx,MAT_POINT,STEP]=...
                Update_ini(STEP,GLOBAL,Disp_field,Mat_state,Int_var,MAT_POINT)
            
            global SOLVER GEOMETRY MATERIAL
            
            BLK=STEP.BLCK;
            %ste=STEP.ste;
            ste_p=STEP.ste_p;
            
            for j=1:GEOMETRY.nodes
                for i=1:GEOMETRY.sp
                     Disp_field.x_a(j,i)=GEOMETRY.x_0(j,i)+...
                         GLOBAL.d((j-1)*GEOMETRY.df+i,ste_p);
                end
             end

            Disp_field.a(:,2)=GLOBAL.a(:,ste_p);
            Disp_field.v(:,2)=GLOBAL.v(:,ste_p);
            Disp_field.d(:,2)=GLOBAL.d(:,ste_p);
            Disp_field.d(:,1)=GLOBAL.d(:,ste_p);

            Int_var.gamma(:,2)=GLOBAL.gamma(:,ste_p);
            Int_var.epsv(:,2)=GLOBAL.epsv(:,ste_p);
            Int_var.Sy(:,2)=GLOBAL.Sy(:,ste_p);
            Int_var.Sy_r(:,2)=GLOBAL.Sy_r(:,ste_p);
            Int_var.eta(:,2)=GLOBAL.eta(:,ste_p);
            Int_var.H(:,2)=GLOBAL.H(:,ste_p);
            Int_var.dgamma(:,2)=GLOBAL.dgamma(:,ste_p);
            
            Mat_state.Sigma(:,2)=GLOBAL.Sigma(:,ste_p);                %STRESS
            
            if SOLVER.SMALL==0
                Mat_state.F(:,1)=GLOBAL.F(:,ste_p);
                Mat_state.F(:,2)=GLOBAL.F(:,ste_p);          %DEFORMATION GRADIENT
                Mat_state.Be(:,2)=GLOBAL.Be(:,ste_p);           %LEFT CAUCHY GREEN
            else
                Mat_state.Es_e(:,2)=GLOBAL.Es(:,ste_p);
                Mat_state.Es(:,2)=GLOBAL.Es_p(:,ste_p)+Mat_state.Es_e(:,2);
                Mat_state.Es(:,1)=Mat_state.Es(:,2);
            end
            
            if SOLVER.UW>0
                Mat_state.pw(:,2)=GLOBAL.pw(:,ste_p);
                if SOLVER.UW==1
                    if SOLVER.SMALL==0
                        Mat_state.Fw(:,2)=GLOBAL.Fw(:,ste_p);
                        Mat_state.Fw(:,1)=GLOBAL.Fw(:,ste_p);
                    else
                        Mat_state.Esw(:,2)=GLOBAL.Esw(:,ste_p);
                    end   
                elseif SOLVER.UW==2
                    Mat_state.dpw(:,2)=GLOBAL.dpw(:,ste_p);
                end
            end
            
            % Initial Stresses   
            MODEL=MATERIAL(BLK).MODEL;
            mmat=MATERIAL(BLK).MAT;
            for i=1:GEOMETRY.mat_points
                mati=GEOMETRY.material(i);
                if MODEL(mati)>=3 && MODEL(mati)<5
                    % P0
                    p0=str2double(mmat{25,mati});
                    ev0=str2double(mmat{23,mati});
                    es0=str2double(mmat{26,mati});
                    if isnan(p0)
                        Int_var.P0(i,1:3)=VECTORS.fill_p0(mmat{25,mati},GLOBAL,i,STEP);
                    else
                        if isnan(ev0)
                            ev0=0;
                        end
                        if isnan(es0)
                            es0=0;
                        end                   
                        Int_var.P0(i,1:3)=[p0 ev0 es0]; 
                    end
                end
            end
            

            %Int_var.P0(:,1)=GLOBAL.Ps(:,1);
           % Mat_state.Sigma(:,3)=GLOBAL.Sigma(:,GLOBAL.final_block(BLK-1));
           Mat_state.Sigma(:,3)=GLOBAL.Sigma(:,1);
            if SOLVER.UW>0
                Mat_state.pw(:,3)=GLOBAL.pw(:,1);
            end
            
            if SOLVER.FRAC
                Mat_state.w(:,1)=GLOBAL.w(:,ste_p);
                Mat_state.status(:,2)=GLOBAL.status(:,ste_p);
                if SOLVER.FRAC>1
                    Mat_state.e_ini(:,1)=GLOBAL.e_ini(:,ste_p);
                end
                STEP.ENERGY.D=GLOBAL.E.D(ste_p,1);
            end
            
            % Recalculate shape functions
            MAT_POINT=shape_function_calculation(0,MAT_POINT,Disp_field);
            
            % Constitutive
            [stiff_mtx,Int_var,Mat_state,STEP]=...
                Constitutive.update(1,STEP,Int_var,Mat_state,MAT_POINT);
            Mat_state.fint(:,2)=Mat_state.fint(:,1);
        end
        
        % Reaction forces
        function [OUTPUT]=reaction(residual,STEP)

            global GEOMETRY SOLVER BOUNDARY
            
            constrains  = BOUNDARY{STEP.BLCK}.constrains;
            
            type=SOLVER.OutputType;
            [number,~]=size(type);
            OUTPUT(1,number)=0;
            
            for i=1:number
                if type(i,1)==1
                    R=0;
                    for j=1:GEOMETRY.df*GEOMETRY.nodes
                        if constrains(j,type(i,2))
                            R=R+residual(j);
                        end
                    end
                    OUTPUT(1,i)=R;
                end
            end

        end
        
        % Plastic Strains in Mat. Points to nodes
        function [Gamma_nds]=Ep2Ep_n(Gamma_tot,MAT_POINT,ste_p)

            global GEOMETRY
            Gamma_nds=zeros(GEOMETRY.nodes,1);

            for i=1:GEOMETRY.mat_points
                sh=MAT_POINT{1}(i).N;
                nds=MAT_POINT{1}(i).near;
                n=length(nds);
                for j=1:n
                    Gamma_nds(nds(j),1)=Gamma_nds(nds(j),1)+sh(j)*Gamma_tot(i,ste_p);
                end
            end

        end
                    
        % Small strain from Def Gradient and Finger tensor
        function [es,es_p]=strains(Mat_state)

            global GEOMETRY SOLVER
            
            dims=GEOMETRY.s_dim;
            [es,es_p]=deal(zeros(GEOMETRY.mat_points*dims,1));
            
            if SOLVER.SMALL==0

                def_G= Mat_state.F;
                b_e  = Mat_state.Be;
                
                dimf=GEOMETRY.f_dim;

                for e=1:GEOMETRY.mat_points
                    [f_v,be]=deal(zeros(dimf,1));
                    for i=1:dimf
                        f_v(i,1)=def_G((e-1)*dimf + i,1);
                        be(i,1)=b_e((e-1)*dimf + i,1);
                    end           
                    [F]=LIB.v2m(f_v);
                    [Be]=LIB.v2m(be);

                    Btot = F*F';

                    [~,R]=chol(Btot);
                    if R==0
                        Etot = logm(Btot)/2;
                        Ee   = logm(Be)/2;
                    else
                        disp('Error in log of B matrx');
                        SOLVER.FAIL=1;
                    end
                    Ep   = Etot-Ee;

                    [ee]=LIB.E2e(Ee);
                    [ep]=LIB.E2e(Ep);
                    for i=1:dims
                        es((e-1)*dims+i,1)=ee(i,1);
                        es_p((e-1)*dims+i,1)=ep(i,1);
                    end
                end
            else
                Es_e= Mat_state.Es_e;
                Es  = Mat_state.Es;
                
                es = Es_e(:,1);
                es_p = Es(:,1) - Es_e(:,1);
            end

        end
        
        % P & Q invariants
        function [P2,Q2]=invar(Ss,e)
            global GEOMETRY
            
            dims=GEOMETRY.s_dim;

            ss=zeros(dims,1);
            for i=1:dims
                ss(dims+1-i)=Ss(e*dims+1-i);
            end

            Sc=LIB.e2E(ss);

            P2=(Sc(1,1)+Sc(2,2)+Sc(3,3))/3;
            s=Sc-P2*eye(3);
            Q2= s(1,1)^2 + s(2,2)^2 + s(3,3)^2 +...
               s(2,1)^2 + s(1,2)^2 + ...
               s(3,1)^2 + s(1,3)^2 + ...
               s(2,3)^2 + s(3,2)^2;
            rj2=Q2*0.5;
            Q2= sqrt(3/2*Q2);
            
                            %Sign
            rj3=0;
            rj3 = rj3 + 3*s(2,1)*s(2,1)*(s(1,1)+s(2,2));
            for i=1:3
                rj3 = rj3 + s(i,i)*s(i,i)*s(i,i);
            end
            rj3=rj3/3;

            rj23=sqrt(rj2)^3;
            if rj23<1.0e-15
                sint3=0;
            else
                sint3 = -3 * sqrt(3) * rj3/2/rj23;
            end
            if sint3<-1
                sint3=-1;
            elseif sint3>1
                sint3=1;
            end
            theta = 1/3*asin(sint3);
            
            if sint3<0
                Q2=-Q2;
            end
        end
        
        % Es & Ev invariants
        function [P2,Q2]=E_invar(Ss,e)
            global GEOMETRY
            
            dims=GEOMETRY.s_dim;

            ss=zeros(dims,1);
            for i=1:dims
                ss(dims+1-i)=Ss(e*dims+1-i);
            end

            Sc=LIB.e2E(ss);

            P2=(Sc(1,1)+Sc(2,2)+Sc(3,3));
            s=Sc-P2*eye(3)/3;
            Q2= s(1,1)^2 + s(2,2)^2 + s(3,3)^2 +...
               s(2,1)^2 + s(1,2)^2 + ...
               s(3,1)^2 + s(1,3)^2 + ...
               s(2,3)^2 + s(3,2)^2;
            rj2=Q2*0.5;
            Q2= sqrt(2/3*Q2);
            P2=-P2;
            
            %Sign
            rj3=0;
            rj3 = rj3 + 3*s(2,1)*s(2,1)*(s(1,1)+s(2,2));
            for i=1:3
                rj3 = rj3 + s(i,i)*s(i,i)*s(i,i);
            end
            rj3=rj3/3;

            rj23=sqrt(rj2)^3;
            if rj23<1.0e-18
                sint3=0;
            else
                sint3 = -3 * sqrt(3) * rj3/2/rj23;
            end
            if sint3<-1
                sint3=-1;
            elseif sint3>1
                sint3=1;
            end
            theta = 1/3*asin(sint3);
            
            if sint3<0
                Q2=-Q2;
            end
        end
        
        
        function [val]=fill_p0(p0,GLOBAL,e,STEP)
            
            global MATERIAL GEOMETRY
            
            Mat=GEOMETRY.material;
            MODEL=MATERIAL(STEP.BLCK).MODEL;
            %MAT=MATERIAL(BLCK).MAT;

            if MODEL(Mat(e))>=3 && MODEL(Mat(e))<4
                disp('Not implemented yet');
                stop;
            elseif MODEL(Mat(e))>=4 && MODEL(Mat(e))<5
                val=VECTORS.p0_PZ(p0,GLOBAL,e,STEP);
            end
            
        end
        
        function [val]=p0_PZ(p0,GLOBAL,e,STEP)
            
            global MATERIAL GEOMETRY
                
            Mat=GEOMETRY.material;
            MAT=MATERIAL(STEP.BLCK).MAT;
            
            if strcmp(p0(1:5),'BLOCK')
                BLCK=str2double(p0(7));
                step=GLOBAL.final_block(BLCK);
            else
                disp('Unrecognized sentence in p0 value');
                stop;
            end
            
            [P,Q]=VECTORS.invar(GLOBAL.Sigma(:,step),e);
            [Ev0,Es1]=VECTORS.E_invar(GLOBAL.Es(:,step),e);
            
            K = MAT{29,Mat(e)};%khar;
            G = MAT{4,Mat(e)};%ghar;
            
            aux=exp(K*Q*Q/6/G/P/P);
            P02=P/aux;
            
            Es2e=P02*Q/P/3/G;
            Es0=Es1-Es2e;
            
            %CHECK
            khar=-K/P02;
            ghar=-G/P02;
            ees=(Es1-Es0);
            eev=0;
            
            p=P02*exp(khar*eev+(3*ghar*khar*(ees^2)/2));
            q=-P02*3*ees*ghar*exp(khar*eev+(3*ghar*khar*(ees^2)/2));
               
            val=[P02,Ev0,Es0];
            
        end
        
    end
 end
    