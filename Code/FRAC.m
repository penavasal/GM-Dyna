 classdef FRAC
    methods(Static)
        
        function [T,A,W,Mat_state,STEP]=process(e,Kt,MAT_POINT,Mat_state,VECS,STEP)

            global MATERIAL GEOMETRY SOLVER

            MODEL=MATERIAL(STEP.BLCK).MODEL;
            Mat=GEOMETRY.material;
            MAT=MATERIAL(STEP.BLCK).MAT;

            Wlist=VECS.Wlist;
            Crit=VECS.Crit;
            Cauchy=VECS.Cauchy;
            T=Cauchy{e};
            W=Wlist(e);
            A=0;

            if SOLVER.FRAC>1
                Max_e=VECS.Max_e;
                FT=VECS.FT;
                %EigenV=VECS.EigenV;
            end
            if Kt==1 || Kt==2 || Kt==4
                Avec=VECS.Avec;
                A=Avec{e};
            end

            vol_e=GEOMETRY.Area(e)*MAT_POINT{1}(e).J;
            h=GEOMETRY.h_ini(e)*sqrt(MAT_POINT{1}(e).J);

            if MODEL(Mat(e),2)==1
                Ceps=MAT{43,Mat(e)}; 
                GC  =MAT{44,Mat(e)};
                if Mat_state.status(e,2)==0 
                    eps_list=MAT_POINT{1}(e).epsilon;

                    if eps_list(1)>0    
                        G=0;
                        % ------ G calculation -----
                        if Crit(e)>0
                            G=Wlist(e)*vol_e;
                            M=vol_e;
                            for i=1:length(eps_list)
                                if Crit(eps_list(i))>0
                                    j=eps_list(i);
                                    vj=GEOMETRY.Area(j)*MAT_POINT{1}(j).J;
                                    G=G+Wlist(eps_list(i))*vj;
                                    M=M+vj;
                                end
                            end
                            G=G*Ceps*h/M;
                        end

                        % ------- G > Gc ?? --------
                        if G > GC
                            Mat_state.status(e,1)=1;
                            T=zeros(3);
                            STEP.ENERGY.D = STEP.ENERGY.D + Wlist(e)*vol_e;
                        end
                    end
                else
                    T=zeros(3);
                end 
            elseif MODEL(Mat(e),2)==2
                %Ceps = MAT{43,Mat(e)};
                wc   = MAT{45,Mat(e)};
                ft   = MAT{46,Mat(e)};
                wk   = MAT{47,Mat(e)};
                fk   = MAT{48,Mat(e)};
                D    = MAT{49,Mat(e)};
                if Mat_state.status(e,2)==0  &&  Mat_state.e_ini(e)==0
                    eps_list=MAT_POINT{1}(e).epsilon;
                    if eps_list(1)>0    
                        f=0;
                        % ------ f calculation -----
                        if Crit(e)>0
                            f=FT(e)*vol_e;
                            M=vol_e;
                            for i=1:length(eps_list)
                                j=eps_list(i);
                                vj=GEOMETRY.Area(j)*MAT_POINT{1}(j).J;
                                if Crit(eps_list(i))>0
                                    f=f+FT(eps_list(i))*vj;
                                    M=M+vj;
                                else
                                    M=M+vj;
                                end
                            end
                            f=f/M;
                        end

                        % ------- G > Gc ?? --------
                        if ft>0 && Crit(e)>0 && f>ft
                            Mat_state.e_ini(e)=Max_e(e);
                        end
                    end
                elseif Mat_state.status(e,2)~=1
                    %chi=(Max_e(e)-E_ini(e))*2*Ceps*h(e)/wc;
                    wn=(Max_e(e)-Mat_state.e_ini(e))*2.5*D;
                    if wn<=wk  
                        chi=wn/wk*(1-fk);
                    else
                        chi=1-fk*(wc-wn)/(wc-wk);
                    end
                    Mat_state.status(e,1)=min(1,max(chi,Mat_state.status(e,2)));
                    T=(1-Mat_state.status(e,1))*T;
                    dd=Wlist(e)*vol_e;
                    eps_list=MAT_POINT{1}(e).epsilon;
                    for i=1:length(eps_list)
                        if Crit(eps_list(i))>0
                            j=eps_list(i);
                            vj=GEOMETRY.Area(j)*MAT_POINT{1}(j).J;
                            dd=dd+Wlist(eps_list(i))*vj;
                        end
                    end

                    STEP.ENERGY.D = STEP.ENERGY.D + ...
                        (Mat_state.status(e,1)-Mat_state.status(e,2))*dd;
                    W=Wlist(e)*(1-Mat_state.status(e,1));
                else
                    T=zeros(3);
                end
            end
        end

        function [MAT_POINT]=eps_nb(MAT_POINT,BLCK)

            global GEOMETRY MATERIAL %Material

            MODEL=MATERIAL(BLCK).MODEL;
            Mat=GEOMETRY.material;
            MAT=MATERIAL(BLCK).MAT;


            h=GEOMETRY.h_ini;

            for e=1:GEOMETRY.mat_points

                Ceps=MAT{43,Mat(e)};

                if MODEL(Mat(e),2)>0

                    C=Ceps*h(e);

                    t=1;
                    for i=1:GEOMETRY.mat_points
                        d=sqrt((MAT_POINT{1}(e).xg(1)-MAT_POINT{1}(i).xg(1))^2+...
                            (MAT_POINT{1}(e).xg(2)-MAT_POINT{1}(i).xg(2))^2);
                        if d<=C & i~=e & MODEL(Mat(i),2)>0 % && Material(e)==Material(i) Different bodies
                            list(t)=i;
                            t=t+1;
                        end
                    end
                    MAT_POINT{1}(e).epsilon=list;

                    if MODEL(Mat(e),2)==3
                        delta95 = MAT{63,Mat(e)};
                        xi95 = MAT{64,Mat(e)};
                        if xi95==0
                            if delta95==0
                                error('Wrong delta and xi 95 properties');
                            else
                                dmax=FRAC.dmax(MAT_POINT,list);
                                xi95=delta95/dmax;
                                MATERIAL(BLCK).MAT(64,M)={xi95};
                            end
                        end
                    end

                    clear list;  

                else
                    MAT_POINT{1}(e).epsilon=0;
                end

            end

        end        


        function [E]=compute_energy(MAT_POINT,STEP,Mat_state,E,v1)


            global GEOMETRY MATERIAL SOLVER 

            step=STEP.ste_p;
            BLCK=STEP.BLCK;

            %MODEL=MATERIAL(BLCK).MODEL;
            mati=GEOMETRY.material;
            %MAT=MATERIAL(BLCK).MAT;


            v2=zeros(GEOMETRY.nodes,1);

            %% Velocity magnitude **********************
            for i=1:GEOMETRY.nodes
                for j=1:GEOMETRY.sp
                    v2(i)=v2(i)+v1(i*GEOMETRY.sp+1-j,1)^2;
                end
                v2(i)=sqrt(v2(i));
            end

            %% Calculation of energy **********************
            for e=1:GEOMETRY.mat_points

                vol_e=GEOMETRY.Area(e)*MAT_POINT{1}(e).J;
                if SOLVER.UW
                    n=1-(1-MATERIAL(BLCK).MAT{16,mati(i)})/MAT_POINT{1}(i).J;
                    dens=n*rho_w+(1-n)*MATERIAL(BLCK).MAT{3,mati(i)};
                else
                    dens=MATERIAL(BLCK).MAT{3,mati(i)}/MAT_POINT{1}(i).J;
                end

                %Strain energy
                if Mat_state.status(e,1)~=0
                    E.W(step,GEOMETRY.body(e))=E.W(step,GEOMETRY.body(e)) + ...
                        Mat_state.w(e) * vol_e;
                end

                %Kinetic energy
                sh = MAT_POINT{1}(e).N;
                nd = MAT_POINT{1}(e).near;
                m  = length(nd);         

                K=0;
                for t1=1:m
                    K = K +dens*vol_e*sh(t1)*v2(nd(t1))*v2(nd(t1))/2;
                end

                E.K(step,GEOMETRY.body(e))=E.K(step,GEOMETRY.body(e)) + K;
            end

    %         %% Calculation of support energy **********************
    %         D=0;
    %         D1=0;
    %         for i=1:nodes
    %             if x_0(i,2)==min(x_0(:,2))
    %                 if (x_0(i,1) <= 72.0*AMP && x_0(i,1) >= 48.0*AMP) || ...
    %                     (x_0(i,1) <= 372.0*AMP && x_0(i,1) >= 348.0*AMP)
    % 
    %                     D=min(D,d(sp*i,ste_p));
    %                     D1=min(D1,d(sp*i,ste_p-1));
    %                 end
    %             end
    %         end
    %         Sup_energy(ste_p)=Sup_energy(ste_p-1)+(React(ste_p)+React(ste_p-1))/2*(D-D1);
        end   

        function [F]=forces(MAT_POINT,STEP,F,v1)


            global GEOMETRY MATERIAL SOLVER 

            dt=STEP.dt;
            step=STEP.ste_p;
            BLCK=STEP.BLCK;

            sp=GEOMETRY.sp;

            %MODEL=MATERIAL(BLCK).MODEL;
            mati=GEOMETRY.material;
            %MAT=MATERIAL(BLCK).MAT;

            vn1=zeros(GEOMETRY.nodes,GEOMETRY.sp);
            vn=zeros(GEOMETRY.nodes,GEOMETRY.sp);

            %% Velocity magnitude **********************
            for i=1:GEOMETRY.nodes
                for j=1:GEOMETRY.sp
                    vn1(i,j)=v1(i*GEOMETRY.sp+1-j,1);
                    vn(i,j) =v1(i*GEOMETRY.sp+1-j,2);
                end
            end

            %% Calculation of energy **********************
            for e=1:GEOMETRY.mat_points

                b=GEOMETRY.body(e);

                vol_e=GEOMETRY.Area(e)*MAT_POINT{1}(e).J;
                if SOLVER.UW
                    n=1-(1-MATERIAL(BLCK).MAT{16,mati(i)})/MAT_POINT{1}(i).J;
                    dens=n*rho_w+(1-n)*MATERIAL(BLCK).MAT{3,mati(i)};
                else
                    dens=MATERIAL(BLCK).MAT{3,mati(i)}/MAT_POINT{1}(i).J;
                end

                sh = MAT_POINT{1}(e).N;
                nd = MAT_POINT{1}(e).near;
                m  = length(nd); 

                % Linear momenta
                L1=0;
                L2=0;
                L1n=0;
                L2n=0;
                for t1=1:m
                    L1 = L1 +dens*vol_e*sh(t1)*vn1(nd(t1),1);
                    L2 = L2 +dens*vol_e*sh(t1)*vn1(nd(t1),2);
                    L1n = L1n +dens*vol_e*sh(t1)*vn(nd(t1),1);
                    L2n = L2n +dens*vol_e*sh(t1)*vn(nd(t1),2);
                end

                % Forces
                F(step,(b-1)*sp+1)=F(step,(b-1)*sp+1) + ...
                    (L1-L1n)/dt;
                F(step,(b-1)*sp+2)=F(step,(b-1)*sp+2) + ...
                    (L2-L2n)/dt;
            end

        end

        function [Edev]=e_deviatoric(Mat_state)

            global GEOMETRY 

            Edev=zeros(GEOMETRY.mat_points,1);

            for e=1:GEOMETRY.mat_points
                dimf=GEOMETRY.f_dim;

                f_v       = zeros(dimf,1);
                f_old     = zeros(dimf,1);
                %Vector F, Fp to Matrixes
                for i=1:dimf
                    f_v(i,1)=Mat_state.F((e-1)*dimf + i,1);
                    f_old(i,1)=Mat_state.F((e-1)*dimf + i,2);
                end           
                [F]=LIB.v2m(f_v);
                [Fk]=LIB.v2m(f_old);

                 B=F*F';
                 Bk=Fk*Fk';
                 E = logm(B)/2;
                 Ek = logm(Bk)/2;
                 dE=E-Ek;

                 Edev(e)=FRAC.invar(dE); 
             end

        end

        function [sy,H,E_ini]=eigendegradation(e,E_dev,E_ini,MAT_POINT,mat95)

            global GEOMETRY

            vol_e=GEOMETRY.Area(e)*MAT_POINT{1}(e).J;
            eps_list=MAT_POINT{1}(e).epsilon;

            dE=E_dev(e)*vol_e;
            M=vol_e;
            for i=1:length(eps_list)
                j=eps_list(i);
                vj=GEOMETRY.Area(j)*MAT_POINT{1}(j).J;

                dE=dE+E_dev(eps_list(i))*vj;
                M=M+vj;
            end
            dE=dE/M;

            E_ini(1)=E_ini(1)+abs(dE);

            if E_ini(2)==0
                sy=mat95(3);
                H=0;
            else
                sy0=mat95(3);
                tau=mat95(1);
                xi95=mat95(2);

                xi=E_ini(1)-E_ini(2);

                sy=tau+(sy0-tau)*exp(-3*xi/xi95);
                H=-3*(sy0-tau)/xi95*exp(-3*xi/xi95);
            end

        end

        function [EP,vec]=Principal(C)
             eps=eig(C);
            [V,~]=eig(C);
            vec=zeros(3);

            [EP(1),i]=max(eps);
            vec(:,1)=V(:,i);
            eps(i)=-1e32;
            [EP(2),i]=max(eps);
            vec(:,2)=V(:,i);
            eps(i)=-1e32;
            [EP(3),i]=max(eps);
            vec(:,3)=V(:,i);
        end

        function [q]=invar(s)
            p=s(1,1)+s(2,2)+s(3,3);

            devs=s-p/3*eye(3);

            [sn]=FRAC.s_j2(devs);

            rj2=sn*sn*0.5;


            rj3=0;
            rj3 = rj3 + 3*devs(2,1)*devs(2,1)*(devs(1,1)+devs(2,2));
            for i=1:3
                rj3 = rj3 + devs(i,i)*devs(i,i)*devs(i,i);
            end
            rj3=rj3/3;

            rj23=sqrt(rj2)^3;
            if 2*rj23<1.0e-18
                sint3=0;
            else
                sint3 = -3 * sqrt(3) * rj3/2/rj23;
            end
            if sint3<-1
                sint3=-1;
            elseif sint3>1
                sint3=1;
            end


            q=sqrt(2/3)*sn;

            if sint3<0
                q=-q;
            end   

        end

        function dmax=dmax(MAT_POINT,list)

            k=0;
            for i=1:length(list)-1
                for j=i+1:length(list)-1
                    k=k+1;
                    d(k)=sqrt((MAT_POINT{1}(j).xg(1)-MAT_POINT{1}(i).xg(1))^2+...
                            (MAT_POINT{1}(j).xg(2)-MAT_POINT{1}(i).xg(2))^2);
                end
            end

            dmax=max(d);
        end
        
        function [q]=s_j2(s)
            %L2-norm of deviatoric stress
            q= s(1,1)^2 + s(2,2)^2 + s(3,3)^2 +...
               s(2,1)^2 + s(1,2)^2 + ...
               s(3,1)^2 + s(1,3)^2 + ...
               s(2,3)^2 + s(3,2)^2;
            q= sqrt(q);
        end

    end
 end