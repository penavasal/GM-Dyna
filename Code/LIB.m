 classdef LIB
    methods(Static)
        
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
            [list]=LIB.S2list(S,field);
            [a,b]=size(list);
            list2=reshape(list,[a*b,1]);
        end
        
        % Convergence and alpha scale (Line search)
        function [CONVER,NORMErec,a,iter]=...
            convergence(r,normr0,NORMErec,toll,iter,imax,a)
        
            global SOLVER
        
            NORMErec(iter,1) = norm(r)/normr0; 
            CONVER=0;
            
            if SOLVER.FAIL==1
                CONVER=1;
            elseif NORMErec(iter,1) < toll  ||  norm(r)< toll
                CONVER=1;
            else
                if iter<imax
                    if iter>imax/2
                        %fprintf('iter %i \n',iter);
                        if std(NORMErec(iter-10:iter-1))<toll/100
                            CONVER=1;
                        end
                    end
                    if (NORMErec(iter)-NORMErec(iter-1))>1e-5 && iter>4
                        f1=NORMErec(iter-1);
                        f2=NORMErec(iter);
                        a=a*a*f1/2/(f2+f1*a-f1);
                        iter = iter-1;
                    else
                        a=1;
                    end
                    if a<1e-9
                        b=min(iter-1,10);
                        if std(NORMErec(iter-b:iter))<toll/10 ||...
                           NORMErec(iter,1) < toll *100
                           CONVER=1;
                        else
                            fprintf('\n No convergence RM \n')
                            stop;
                        end
                    end
                    
                elseif iter == imax

                    if NORMErec(iter,1)<toll*10
                        CONVER=1;
                    %elseif std(normr0*NORMErec(iter-10:iter-1))<toll/10
                    elseif std(NORMErec(iter-10:iter-1))<toll/2
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