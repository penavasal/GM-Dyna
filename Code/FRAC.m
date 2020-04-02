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
            Ceps=MAT(43,Mat(e)); 
            GC  =MAT(44,Mat(e));
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
            %Ceps = MAT(43,Mat(e));
            wc   = MAT(45,Mat(e));
            ft   = MAT(46,Mat(e));
            wk   = MAT(47,Mat(e));
            fk   = MAT(48,Mat(e));
            D    = MAT(49,Mat(e));
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

            Ceps=MAT(43,Mat(e));

            if MODEL(Mat(e),2)>0

                C=Ceps*h(e);

                t=1;
                for i=1:GEOMETRY.mat_points
                    d=sqrt((MAT_POINT{1}(e).xg(1)-MAT_POINT{1}(i).xg(1))^2+...
                        (MAT_POINT{1}(e).xg(2)-MAT_POINT{1}(i).xg(2))^2);
                    if d<=C && i~=e % && Material(e)==Material(i) Different bodies
                        list(t)=i;
                        t=t+1;
                    end
                end
                MAT_POINT{1}(e).epsilon=list;
                clear list;  

            else
                MAT_POINT{1}(e).epsilon=0;
            end

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
                n=1-(1-MATERIAL(BLCK).MAT(16,mati(i)))/MAT_POINT{1}(i).J;
                dens=n*rho_w+(1-n)*MATERIAL(BLCK).MAT(3,mati(i));
            else
                dens=MATERIAL(BLCK).MAT(3,mati(i))/MAT_POINT{1}(i).J;
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
                n=1-(1-MATERIAL(BLCK).MAT(16,mati(i)))/MAT_POINT{1}(i).J;
                dens=n*rho_w+(1-n)*MATERIAL(BLCK).MAT(3,mati(i));
            else
                dens=MATERIAL(BLCK).MAT(3,mati(i))/MAT_POINT{1}(i).J;
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
    
    
    end
 end