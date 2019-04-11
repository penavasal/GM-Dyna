

 classdef AUX
    methods(Static)
        
        % P & Q invariants
        function [P2,Q2]=invar(Ss,e)

            ss=zeros(4,1);
            for i=1:4
                ss(5-i)=Ss(e*4+1-i);
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

            E=zeros(3,3);

            %Build matrix

            E(1,1)=e(1);
            E(2,2)=e(2);
            E(1,2)=e(4);
            E(2,1)=e(4);
            E(3,3)=e(3);

        end
        
        % STRESS - STRAIN tensor to vector
        function [e]=E2e(E)
            e=zeros(4,1);

            %Build vector
            e(1)=E(1,1);
            e(2)=E(2,2);
            e(3)=E(3,3);
            e(4)=E(1,2);%+E(2,1);
        end
        
        % DEF GRADIENT tensor to vector
        function [f]=m2v(F,sp)

            f=zeros(sp*sp+1,1);

            %Build vector
            for i=1:sp
                for j=1:sp
                   f((i-1)*sp+j,1)=F(i,j);
                end
            end
            f(sp*sp+1,1)=F(i+1,j+1);
        end
        
        % DEF GRADIENT vector to tensor
        function [F]=v2m(def_G,sp)
   
            F=zeros(sp+1,sp+1);

            %Build matrix
            for i=1:sp
                for j=1:sp
                    F(i,j)=def_G((i-1)*sp+j,1);
                end
            end
            F(i+1,j+1)=def_G(sp*sp+1,1);
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

            [es,es_p]=deal(zeros(GEOMETRY.mat_points*4,1));
            [f_v,be]=deal(zeros(5,1));

            for e=1:GEOMETRY.mat_points
                for i=1:5
                    f_v(i,1)=def_G((e-1)*5 + i,1);
                    be(i,1)=b_e((e-1)*5 + i,1);
                end           
                [F]=AUX.v2m(f_v,GEOMETRY.sp);
                [Be]=AUX.v2m(be,GEOMETRY.sp);

                Btot = F*F';
                Etot = logm(Btot)/2;
                Ee   = logm(Be)/2;
                Ep   = Etot-Ee;

                [ee]=AUX.E2e(Ee);
                [ep]=AUX.E2e(Ep);
                for i=1:4
                    es((e-1)*4+i,1)=ee(i,1);
                    es_p((e-1)*4+i,1)=ep(i,1);
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
        
        % Structure to list
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
    end
 end
 
 
 
 
 
 