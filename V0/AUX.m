

 classdef AUX
    methods(Static)
        
        % P & Q invariants
        function [P2,Q2]=invar(Ss,e)
            global GEOMETRY
            
            dims=GEOMETRY.s_dim;

            ss=zeros(dims,1);
            for i=1:dims
                ss(dims+1-i)=Ss(e*dims+1-i);
            end

            Sc=AUX.e2E(ss);

            P2=(Sc(1,1)+Sc(2,2)+Sc(3,3))/3;
            s=Sc-P2*eye(3);
            Q2= s(1,1)^2 + s(2,2)^2 + s(3,3)^2 +...
               s(2,1)^2 + s(1,2)^2 + ...
               s(3,1)^2 + s(1,3)^2 + ...
               s(2,3)^2 + s(3,2)^2;
            Q2= sqrt(3/2*Q2);
        end
        
        % STRESS - STRAIN vector to tensor
        function [E]=e2E(e)

            global GEOMETRY 
            
            if GEOMETRY.sp==1
                E=e(1);
            else
                E=zeros(3,3);
                %Build matrix
                if GEOMETRY.sp==2
                    E(1,1)=e(1);
                    E(2,2)=e(2);
                    E(1,2)=e(4);
                    E(2,1)=e(4);
                    E(3,3)=e(3);
                else
                    E(1,1)=e(1);
                    E(2,2)=e(2);
                    E(1,2)=e(4);
                    E(2,1)=e(4);
                    E(1,3)=e(5);
                    E(3,1)=e(5);
                    E(2,3)=e(6);
                    E(3,2)=e(6);
                    E(3,3)=e(3);
                end
            end
        end
        
        % STRESS - STRAIN tensor to vector
        function [e]=E2e(E)
            
            global GEOMETRY
            
            e=zeros(GEOMETRY.s_dim,1);
            
            %Build vector
            e(1)=E(1,1);
            if GEOMETRY.sp>1
                e(2)=E(2,2);
                e(3)=E(3,3);
                e(4)=E(1,2);
                if GEOMETRY.sp>2
                	e(5)=E(1,3);
                    e(6)=E(2,3);
                end
            end              
        end
        
        % STRESS - STRAIN vector to tensor
        function [E]=e2E_in(e)

            global GEOMETRY 
            
            if GEOMETRY.sp==1
                E=e(1);
            else
                E=zeros(3,3);
                %Build matrix
                if GEOMETRY.sp==2
                    E(1,1)=e(1);
                    E(2,2)=e(2);
                    E(1,2)=e(3);
                    E(2,1)=e(3);
                    E(3,3)=e(4);
                else
                    E(1,1)=e(1);
                    E(2,2)=e(2);
                    E(1,2)=e(3);
                    E(2,1)=e(3);
                    E(1,3)=e(5);
                    E(3,1)=e(5);
                    E(2,3)=e(6);
                    E(3,2)=e(6);
                    E(3,3)=e(4);
                end
            end
        end
        
        % STRESS - STRAIN tensor to vector
        function [e]=E2e_in(E)
            
            global GEOMETRY
            
            e=zeros(GEOMETRY.s_dim,1);
            
            %Build vector
            e(1)=E(1,1);
            if GEOMETRY.sp>1
                e(2)=E(2,2);
                e(3)=E(1,2);
                e(4)=E(3,3);
                if GEOMETRY.sp>2
                	e(5)=E(1,3);
                    e(6)=E(2,3);
                end
            end
                    
        end
        
        % DEF GRADIENT tensor to vector
        function [f]=m2v(F)

            global GEOMETRY
            
            sp=GEOMETRY.sp;
            dimf=GEOMETRY.f_dim;
            
            f=zeros(dimf,1);

            %Build vector
            for i=1:sp
                for j=1:sp
                   f((i-1)*sp+j,1)=F(i,j);
                end
            end
            if sp==2
                f(sp*sp+1,1)=F(i+1,j+1);
            end
        end
        
        % DEF GRADIENT vector to tensor
        function [F]=v2m(def_G)
   
            global GEOMETRY
            
            sp=GEOMETRY.sp;
            dimf=GEOMETRY.f_dim;
            
            if sp==1
                F=zeros(sp,sp);
            else
                F=zeros(3,3);
            end

            %Build matrix
            for i=1:sp
                for j=1:sp
                    F(i,j)=def_G((i-1)*sp+j,1);
                end
            end
            if sp==2
                F(3,3)=def_G(dimf,1);
            end
        end
        
        % Reaction forces
        function [OUTPUT]=reaction(residual,OUTPUT)

            global GEOMETRY

            for i=1:OUTPUT.number
                if OUTPUT.type(i,1)==2
                    R=0;
                    for j=1:GEOMETRY.df*GEOMETRY.nodes
                        if OUTPUT.ref_list(j,i)==1
                            R=R+residual(j);
                        end
                    end
                    OUTPUT.inst(i)=R;
                end
            end

        end
        
        % Small strain from Def Gradient and Finger tensor
        function [es,es_p]=strains(def_G,b_e)

            global GEOMETRY

            dimf=GEOMETRY.f_dim;
            dims=GEOMETRY.s_dim;
            
            [es,es_p]=deal(zeros(GEOMETRY.mat_points*dims,1));
            [f_v,be]=deal(zeros(dimf,1));

            for e=1:GEOMETRY.mat_points
                for i=1:dimf
                    f_v(i,1)=def_G((e-1)*dimf + i,1);
                    be(i,1)=b_e((e-1)*dimf + i,1);
                end           
                [F]=AUX.v2m(f_v);
                [Be]=AUX.v2m(be);

                Btot = F*F';
                
                [~,R]=chol(Btot);
                if R==0
                    Etot = logm(Btot)/2;
                    Ee   = logm(Be)/2;
                else
                    disp('Error in log of B matrx');
                    pause
                end
                Ep   = Etot-Ee;

                [ee]=AUX.E2e(Ee);
                [ep]=AUX.E2e(Ep);
                for i=1:dims
                    es((e-1)*dims+i,1)=ee(i,1);
                    es_p((e-1)*dims+i,1)=ep(i,1);
                end
            end

        end
        
        % Plastic Strains in Mat. Points to nodes
        function [Gamma_nds]=Ep2Ep_n(Gamma_tot,MAT_POINT,ste_p)

            global GEOMETRY
            Gamma_nds=zeros(GEOMETRY.nodes,1);

            for i=1:GEOMETRY.mat_points
                sh=MAT_POINT(i).N;
                nds=MAT_POINT(i).near;
                n=length(nds);
                for j=1:n
                    Gamma_nds(nds(j),1)=Gamma_nds(nds(j),1)+sh(j)*Gamma_tot(i,ste_p);
                end
            end

        end
        
        % Structure to list
        function [list]=S2list(S,field)
            [~,b]=size(S);
            %list=zeros(b,1);
            for i=1:b
                list(i,:) = getfield(S(i),field);
            end
        end
        
        % list to Structure
        function [S]=list2S(S,field,list)
            [a,~]=size(list);
            %list=zeros(b,1);
            for i=1:a
                S(i) = setfield(S(i),field,list(i,:));
            end
        end
        
        % In or out search function
        function [I]=IoO(X,x_a,list)
    
            [~,nCC]=size(list);

            d=zeros(nCC,1);
            I=1;
            for i=1:nCC
                x0=x_a(list(i),:);
                if i==nCC
                    xf=x_a(list(1,1),:);
                else
                    xf=x_a(list(1,i+1),:);
                end
                M=(xf(2)-x0(2))/(xf(1)-x0(1));
                if isinf(M) || abs(M)>1e2
                    if (xf(2)-x0(2))>0
                        M=1;
                    else
                        M=-1;
                    end
                    d(i)=M*(X(1,1)-x0(1,1));
                else
                    if (xf(1)-x0(1))>0
                        M0=1;
                    else
                        M0=-1;
                    end
                    rM=sqrt(M^2+1);
                    d(i)=M0*(M*X(1,1)-X(1,2)+x0(1,2)-M*x0(1,1))/rM;
                end

                if i~=1 && d(1)*d(i)<=0
                    I=0;
                end
            end
        end
        
        % Reshape structure to save
        function [list2]=reshape_S2list(S,field)
            [list]=AUX.S2list(S,field);
            [a,b]=size(list);
            list2=reshape(list,[a*b,1]);
        end
        
        % Convergence and alpha scale (Line search)
        function [CONVER,NORMErec,a,iter]=...
            convergence(r,normr0,NORMErec,toll,iter,imax,a)
        
            NORMErec(iter,1) = norm(r)/normr0; 
            CONVER=0;
            if NORMErec(iter,1) < toll  ||  norm(r)< toll
                CONVER=1;
            else
                if iter<imax
                    if iter>imax/2
                        fprintf('iter %i \n',iter);
                        if std(NORMErec(iter-10:iter-1))<toll/100
                            CONVER=1;
                        end
                    else
                    
                        if NORMErec(iter)>NORMErec(iter-1) && iter>4
                            f1=NORMErec(iter-1);
                            f2=NORMErec(iter);
                            a=a*a*f1/2/(f2+f1*a-f1);
                            iter = iter-1;
                        else
                            a=1;
                        end
                    end
                    if a<1e-9
                        fprintf('\n No convergence RM \n')
                        stop;
                    end
                    
                elseif iter == imax

                    if NORMErec(iter,1)<toll*10
                        CONVER=1;
                    %elseif std(normr0*NORMErec(iter-10:iter-1))<toll/10
                    elseif std(NORMErec(iter-10:iter-1))<toll/10
                        CONVER=1;
                    else
                        fprintf('\n No convergence RM \n')
                        stop;
                    end
                end         
            end
            
        end

        
    end
 end
 
 
 
 
 
 