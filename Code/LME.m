
 classdef LME
    methods(Static)
        
        % File: LME
        %   
        % Date:
        %   Version 2.0   10.04.2019

        function [MAT_POINT]=near(i,x_a,x,h,MAT_POINT,sh,ph)

            global GEOMETRY LME_param 

            e=MAT_POINT(i).element; %Element where is the material point
            
            no=MAT_POINT(i).near_p;
            nn=length(no);
            if no==0
                nn=0;
            end

            ndim=GEOMETRY.sp;

            target_zero= LME_param(sh,2);

            nbg=LME_param(sh,8);
            % Range of search
            range_aux=h(i)*sqrt(-log(target_zero)/MAT_POINT(i).w);
            range=max(range_aux, nbg*h(i));
                       
            % Nodes where seeking
            nds_search = MAT_POINT(i).near;
            if nbg==1
                el_near=GEOMETRY.element_near{e};
                for j=1:length(el_near)
                    if ph==2 && (strcmp(GEOMETRY.ELEMENT,'Q8P4') || ...
                        strcmp(GEOMETRY.ELEMENT,'Q8P4-4'))
                        elem=GEOMETRY.elem_c(el_near(j),:);
                    else
                        elem=GEOMETRY.elem(el_near(j),:);
                    end
                    C = setdiff(elem,nds_search);
                    if C
                        nds_search=cat(2,nds_search,C);
                    end
                end
            elseif nbg==2
                el_near=GEOMETRY.element_near{e};
                for k=1:length(el_near)
                    el_near_2=GEOMETRY.element_near{el_near(k)};
                    for j=1:length(el_near_2)
                        if ph==2 && (strcmp(GEOMETRY.ELEMENT,'Q8P4') || ...
                            strcmp(GEOMETRY.ELEMENT,'Q8P4-4'))
                            elem=GEOMETRY.elem_c(el_near_2(j),:);
                        else
                            elem=GEOMETRY.elem(el_near_2(j),:);
                        end
                        C = setdiff(elem,nds_search);
                        if C
                            nds_search=cat(2,nds_search,C);
                        end
                    end
                end
            elseif nbg==0
                nds_search=[];
            else
                disp('Error, not defined so many neighborhood ranges')
                stop
            end
                        
            % Coordinates of the possible nodes
            x_find=zeros(length(nds_search),ndim);
            for j=1:length(nds_search)
                x_find(j,:)=x_a(nds_search(j),:);
            end
            
            %Naive find in range
            dist=zeros(length(nds_search),1);
            for id=1:ndim
                dist=dist + (x_find(:,id)-x(id)).^2;
            end
            dist=sqrt(dist);

            if ph==2 && (strcmp(GEOMETRY.ELEMENT,'Q8P4') || ...
                    strcmp(GEOMETRY.ELEMENT,'Q8P4-4'))
                corner=GEOMETRY.elem_c;
                near_=GEOMETRY.elem_c(e,:);
            else
                corner=GEOMETRY.elem;
                near_=GEOMETRY.elem(e,:);
            end
            T=0;
            for j=1:length(nds_search)
              %if dist(j)<range(i)
              if  (dist(j)<range) && not(ismember(nds_search(j),near_))...
                      && any(ismember(corner(:),nds_search(j)))%...
                      % &&(Mat(i)==Mat_nds(nds_search(j)))
                  t=0;
                  for k=1:nn
                    if nds_search(j)==no(k)
                        t=1;
                        T=1;
                    end
                  end   
                  
                  if t==0
                    near_=cat(2,nds_search(j),near_);
                  end
              end
            end
            MAT_POINT(i).near=near_;
            if T==1
                plot_nb(i,near_,x_a,GEOMETRY.elem_c,1,0);
            end
        end

        function [MAT_POINT,FAIL]=shapef(i,x_a,h,x,MAT_POINT,sh)
        
    % Calculation of the local maximum-entropy shape functions with Newton's
    % method as described in section 4.2 of [1]
    %
    % INPUT
    % =====
    % ndim:     spacial dimension
    % x_a:      nodal set (nnode, ndim)
    % beta_:    value of the thermalization parameter at each sample point
    % x_sample: sample point where the shape functions are evaluated
    % near_:    near_{i} contains the list of nodes with non-zero shape funcion
    %           at the i-th sample point
    % TolLag:   tolerance for Newton's iterations
    %
    % OUTPUT
    % ======
    % p:        p{i} contains the values of the shape functions corresponding
    %           to the nodes in near_{i} at the i-th sample point
    % dp:       dp{i} contains the values of the ndim spacial derivatives of 
    %           the shape functions corresponding to the nodes in near_{i} at 
    %           the i-th sample point
    %
    % Reference:
    % [1] Marino Arroyo and Michael Ortiz, "Local maximum-entropy approximation
    %     schemes: a seamless bridge between finite elements and meshfree methods", 
    %     International Journal for Numerical Methods in Engineering, 65:2167-2202 (2006). 

            global LME_param GEOMETRY
            
            ndim = GEOMETRY.sp;
    
            TolLag= LME_param(sh,3);
            Nelder=LME_param(sh,4); 
    
            FAIL=0;

            gamma_=MAT_POINT(i).w;
            beta=gamma_/(h(i))^2;
            near= MAT_POINT(i).near;
            lam = MAT_POINT(i).xi;

            %lam=lama';

            if Nelder
              [p_a,J,cnt,lam,tol]=LME.Nelder_mead(x,x_a,beta,near,TolLag,lam,ndim);
            else
              [p_a,J,cnt]=LME.Newton_raphson(x,x_a,beta,near,TolLag,lam,ndim);
            end

            if cnt>=250
              fprintf('Convergence problem in element %i with tol %i \n'...
                  ,i,tol);
            end

            MAT_POINT(i).xi=lam;
            MAT_POINT(i).N=p_a; 

            dp=zeros(length(near),ndim);
            if isnan(p_a)
              fprintf('Not a number in element %i \n',i);
              FAIL=1;
            elseif isnan(J)
              fprintf('Not a number in element %i \n',i);
              FAIL=1;
            else
              %Spacial Gradients
              for ia=1:length(near)
                if ndim==1
                  dp(ia,1)= -(p_a(ia) * (x-x_a(near(ia))') )/J(1,1);      
                else
                  dp(ia,:)= -J\(p_a(ia) * (x(:)-x_a(near(ia),:)') );
                end
              end
            end
            MAT_POINT(i).B=dp;
        end

        function [lambda]=first_lambda(i,elem,x_sample,x_a,beta)

            [~,bg]=size(elem);
            [~,sp]=size(x_a);

            x=x_sample(i,:);

            % Initialize lambda
            y=zeros(bg,sp);
            for j=1:bg
                for k=1:sp
                    y(j,k)=x(k)-x_a(elem(j),k);
                end
            end
            A(1,1)=y(2,1)-y(1,1);
            A(1,2)=y(2,2)-y(1,2);
            A(2,1)=y(3,1)-y(1,1);
            A(2,2)=y(3,2)-y(1,2);
            B(1,1)=-beta*(y(1,1)^2+y(1,2)^2-(y(2,1)^2+y(2,2)^2));
            B(2,1)=-beta*(y(1,1)^2+y(1,2)^2-(y(3,1)^2+y(3,2)^2));
            lam=A\B;

            lambda=lam';
            
%             N0=[0.3333333; 0.3333333];
%             L=[0 0];
%             R=N0;
%             rtol=1e-4;
%             while (norm(R)>rtol)
%                 [~,~,~,N]=LME.Gamma_(sp,x_a,x,beta,L,elem);
%                 R=N(1:2)-N0;
%                 J=[N(1)*y(1,1) N(1)*y(1,2);
%                    N(2)*y(2,1) N(2)*y(2,2)];
%                 dL=J\R;
%                 L=L-dL';
%             end
%             
%             lambda=L';

        end

        function [p_a,J,niter]=Newton_raphson(x,x_a,beta,near,TolLag,lam,ndim)


          R=10*ones(1,ndim);
          niter=0;
          I=eye(ndim);

          %Newton iteration
          while (norm(R)>TolLag)
            [gam,R,J,p_a]=LME.Gamma_(ndim,x_a,x,beta,lam,near);
            if (abs(rcond(J)))<1e-8
              iflag=1;
              disp('Newton Failed, near to singular matrix')
            end
            dlam=-J\R;
            lam=lam+dlam';
            niter=niter+1;
            if (niter>100) 
              i, niter
              disp('Newton Failed 2, no convergence in 100 iterations')
            end
          end

        end

        function [p_a,J,cnt,lam,TOL]=Nelder_mead(x,x_a,beta,near,TolLag,lam,ndim)
            %
            max1=250;
            min1=5;
            %
            rho=1;
            xi=2;
            gam=0.5;
            sig=0.5;
            epsilon=1e-2;

            pts=ndim+1;
            f(pts)=0;

            %Lambda


            L(1,:)=[lam(1),lam(2)/10];
            L(2,:)=[lam(1),lam(2)];
            L(3,:)=[lam(1)/10,lam(2)];

            for i=1:pts
                lam=L(i,:);
                [f(i),~,~,~]=LME.Gamma_(ndim,x_a,x,beta,lam,near);
            end


            % 1- Ordenar

            [LL,f_m,f_mm,f_w]=LME.order(L,f);

            %%%%% ROTATION ?? %%%%%
            ROT=1;
            iter=1;
            while ROT==1 && iter<100
                [Vn,diam]=LME.von(LL);
                if Vn < epsilon
                    ROT=1;
                    [LL,f_m,f_mm,f_w]=LME.rota2D(LL,near,x_a,x,beta,diam);
                else
                    ROT=0;
                end
                iter=iter+1;
            end
            %%%%%%%%%%%%%%%

            cnt=1;
            tol(cnt)=10;
            while (tol(cnt)>TolLag && cnt<max1) || cnt<min1

                %2 - Reflexion
                sum=zeros(1,ndim);
                for i=1:ndim
                    sum=sum+LL(i,:);
                end
                x_bar=sum/ndim;
                x_r=x_bar+rho*(x_bar-LL(3,:));

                [f_r,~,~,~]=LME.Gamma_(ndim,x_a,x,beta,x_r,near);

                if f_r<f_m   
                    %3-Expansion
                    x_e=x_bar+xi*(x_r-x_bar);
                    %x_e=x_r+xi*(x_r-x_bar);
                    [f_e,~,~,~]=LME.Gamma_(ndim,x_a,x,beta,x_e,near);

                    if f_e<f_m
                        LL(3,:)=LL(2,:);
                        LL(2,:)=LL(1,:);
                        LL(1,:)=x_e;
                        f_w=f_mm;
                        f_mm=f_m;
                        f_m=f_e;
                    else
                        LL(3,:)=LL(2,:);
                        LL(2,:)=LL(1,:);
                        LL(1,:)=x_r;
                        f_w=f_mm;
                        f_mm=f_m;
                        f_m=f_r; 
                    end

                elseif f_r<f_mm
                    LL(3,:)=LL(2,:);
                    LL(2,:)=x_r;
                    f_w=f_mm;
                    f_mm=f_r;

                else
                    cont=0;
                    if f_r<f_w
                        %4- Contraccion fuera
                        x_c=x_bar+gam*(x_r-x_bar);
                        [f_c,~,~,~]=LME.Gamma_(ndim,x_a,x,beta,x_c,near);
                        if f_c<f_r
                            cont=1;
                        end

                    else
                        %4- Contraccion dentro
                        x_c=x_bar-gam*(x_bar-LL(3,:));
                        [f_c,~,~,~]=LME.Gamma_(ndim,x_a,x,beta,x_c,near);
                        if f_c<f_w
                            cont=1;
                        end
                    end

                    if cont
                        if f_c<f_m
                            LL(3,:)=LL(2,:);
                            LL(2,:)=LL(1,:);
                            LL(1,:)=x_c;
                            f_w=f_mm;
                            f_mm=f_m;
                            f_m=f_c;

                        elseif f_c<f_mm
                            LL(3,:)=LL(2,:);
                            LL(2,:)=x_c;
                            f_w=f_mm;
                            f_mm=f_c;
                        else
                            LL(3,:)=x_c;
                            f_w=f_c;
                        end

                    else
                        %5- Encogimiento

                        for i=1:pts
                            LL(i,:)=LL(1,:)+sig*(LL(i,:)-LL(1,:));
                            lam=LL(i,:);
                            [f(i),~,~,~]=LME.Gamma_(ndim,x_a,x,beta,lam,near);
                        end
                        f(1)=f_m;
                        [LL,f_m,f_mm,f_w]=LME.order(LL,f);
                    end
                end
                cnt=cnt+1;
                %%%%% ROTATION ?? %%%%%
                ROT=1;
                iter=1;
                while ROT==1 && iter<100
                    [Vn,diam]=LME.von(LL);
                    if Vn < epsilon
                        ROT=1;
                        [LL,f_m,f_mm,f_w]=LME.rota2D(LL,near,x_a,x,beta,diam);
                    else
                        ROT=0;
                    end
                    iter=iter+1;
                end
                %%%%%%%%%%%%%%%
                tol(cnt)=abs(f_w-f_m);
                if cnt>=230
                    cnt;
                end   
            end    
            TOL=tol(cnt);
            lam=LL(1,:);
            [~,~,J,p_a]=LME.Gamma_(ndim,x_a,x,beta,lam,near);
        end

        function [gam,dgam,hgam,p_a]=Gamma_(ndim,x_a,x,beta,lam,near)

            if ndim==1
                temp=exp(-beta*((x-x_a(near)).^2) ...
                    + lam*(x-x_a(near)) );
                Z=sum(temp);
                p_a=temp/Z;
                gam=log(Z);
                temp1(1:length(near))=(x-x_a(near)).*p_a;
                dgam=sum(temp1(:));
                hgam = sum( p_a.*(x-x_a(near)).^2   ) - dgam*dgam;
            else
                sum1=0;
                sum2=0;
                for id=1:ndim
                    sum1=sum1 + (x(id)-x_a(near,id)).^2;
                    sum2=sum2 + lam(id)*(x(id)-x_a(near,id));
                end

                temp=exp(-beta*sum1 + sum2);
                Z=sum(temp);
                p_a=temp/Z;
                gam=log(Z);

                dgam=zeros(ndim,1);
                hgam=zeros(ndim,ndim);
                for id=1:ndim
                    dgam(id)=sum( (x(id)-x_a(near,id)).*p_a(:) );
                end
                for id=1:ndim
                    for jd=1:ndim
                        hgam(id,jd)=sum( p_a(:).* ...
                            (x(id)-x_a(near,id)).*(x(jd)-x_a(near,jd)) )  ...
                            - dgam(id)*dgam(jd);
                    end
                end
            end

        end

        function [LL,f_m,f_mm,f_w]=order(L,f)
            [pts,~]=size(L);

            [f_m, lo] = min(f);            
            [f_w, hi] = max(f);

            for i=1:pts
                if i~=lo && i~=hi
                    mm=i;
                end
            end
            f_mm=f(mm);
            LL(1,:)=L(lo,:);
            LL(2,:)=L(mm,:);
            LL(3,:)=L(hi,:);
        end

        function [Vn,diam]=von(LL)
            [n,sp]=size(LL);
            mat=ones(n);

            mat(:,1:2)=LL;
            DT=abs(det(mat)/factorial(sp));

            dist(n,n)=0;
            for i=1:n
                for k=1:n
                    for j=1:sp
                        dist(i,k)=dist(i,k)+(LL(i,j)-LL(k,j))^2;
                    end
                    dist(i,k)=sqrt(dist(i,k));
                end
            end
            diam=max(max(dist));

            Vn=DT/diam^sp;
        end

        function [LL2,f_m,f_mm,f_w]=rota2D(LL,near,x_a,x,beta,diam)

            [~,sp]=size(LL);

            DIS=max(diam/50,eps*2);


            LL1(1,:)=LL(1,:);
            LL2(1,:)=LL(1,:);

            LL1(2,1)=LL(2,1)+DIS;
            LL1(2,2)=LL(2,2);
            LL1(3,1)=LL(3,1);
            LL1(3,2)=LL(3,2)-DIS;

            LL2(2,1)=LL(2,1)-DIS;
            LL2(2,2)=LL(2,2);
            LL2(3,1)=LL(3,1);
            LL2(3,2)=LL(3,2)+DIS;

            f_1(1)=1e10;
            f_2(1)=1e10;
            [f_1(2),~,~,~]=LME.Gamma_(sp,x_a,x,beta,LL1(2,:),near);
            [f_1(3),~,~,~]=LME.Gamma_(sp,x_a,x,beta,LL1(3,:),near);
            [f_2(2),~,~,~]=LME.Gamma_(sp,x_a,x,beta,LL2(2,:),near);
            [f_2(3),~,~,~]=LME.Gamma_(sp,x_a,x,beta,LL2(3,:),near);

            f1=min(f_1);
            f2=min(f_2);

            if f1<f2
                LL2=LL1;
                f_2=f_1;
            end
            [f_2(1),~,~,~]=LME.Gamma_(sp,x_a,x,beta,LL1(1,:),near);

            [LL2,f_m,f_mm,f_w]=LME.order(LL2,f_2);
        end

        function SEP=read

                global LME_param SOLVER


                % As described in [1], gamma measures the degree of locality, the larger,
                % the closer to Delaunay finite elements

                gamma_lme=3.5;
                gamma_top=0.9;

                % This sets the numerical threshold for the support of the shape functions
                target_zero=1.e-6; 
                % This is the tolerance for the Newton iteration to compute the shape
                % functions
                TolLag=max(2*eps,1.e-18);
                % Nelder Mead method to obtain minimum lambda in shape functions
                Nelder=1;

                %Tolerance to make new re-mapping or not and proportion to do it
                tol_search = 0.5;       % Optimum [0.4-0.7]
                prop = 0.1;             % Proportion of gamma decreasing


                str= SOLVER.TYPE{3};
                
                
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % Read LME.txt
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fid = fopen(str, 'rt'); % opción rt para abrir en modo texto
                formato = '%s %s %s %s %s %s %s %s %s'; % formato de cada línea 
                data = textscan(fid, formato, 'HeaderLines', 1);
                a = data{1};
                % Convertir a vector numérico
                b1= cellfun(@str2num, data{2}, 'UniformOutput', false);
                b2= data{2};
                b3= data{3};
                [l,~] = size(b1);
                b=zeros(l,1);
                for i=1:l
                    if b1{i}
                        b(i)=b1{i};
                    else
                        b(i)=0;
                    end
                end
                t=0;
                
                shfs=str2double(b2{1});
                if isnan(shfs)
                    shfs=0;
                end
                while shfs==0
                    t=t+1;
                    if strcmp(a{t},'SHAPE_FUNCTIONS')
                        shfs=str2double(b2{t});
                    end
                end
                [phases,~]=size(SOLVER.PHASES);
                
                SEP(shfs).l=0;
                SEP(shfs).h=0;
                
                M=0;
                while (M<shfs+1) &&(t<l)
                    t=t+1;
                    s1=a{t};
                    switch s1
                        case '//'
                            continue
                        case 'PHASE'
                            if (M+1)>shfs
                                break;
                            elseif M>0
                                LME_param(M,1)=gamma_lme;
                                LME_param(M,2)=target_zero;
                                LME_param(M,3)=max(2*eps,TolLag);
                                LME_param(M,4)=Nelder;
                                LME_param(M,5)=tol_search;
                                LME_param(M,6)=prop;
                                LME_param(M,7)=gamma_top;
                                LME_param(M,8)=nb_grade;
                            end
                            ss=0;
                            M=M+1;
                            ph=b2{t};
                            
                            for i=1:length(ph)
                                for k=1:phases
                                    if strcmp(ph(i),SOLVER.PHASES{k,1})
                                        SOLVER.PHASES{k,2}=M;
                                    end
                                end
                            end
                            continue
                        case 'GAMMA_LME'
                            gamma_lme=b(t); 
                            continue
                        case 'GAMMA_TOP'
                            gamma_top=b(t);
                            continue
                        case 'TARGET_ZERO'
                            target_zero=b(t);
                            continue
                        case 'TOL_LAG'
                            TolLag=max(b(t),TolLag);  
                            continue
                        case 'TOL_SEARCH'
                            tol_search=b(t);  
                            continue
                        case 'WRAPPER'
                            if strcmp(b2{t},'NELDER') || strcmp(b2{t},'NELDER_MEAD')
                                Nelder=1;
                            elseif strcmp(b2{t},'NR') || strcmp(b2{t},'NEWTON_RAPHSON')
                                Nelder=0;
                            end
                            continue
                        case 'PROPORTION'
                            prop=b(t);  
                            continue
                        case 'NEIGHBORHOOD_GRADE'
                            nb_grade=b(t);  
                            continue
                        case 'SEPARATION'
                            ss=ss+1;
                            SEP(M).l(ss)=b(t);
                            SEP(M).h(ss)=str2double(b3{t});
                            continue
                        otherwise
                            fprintf('Error, unrecognized parameter: %s !!\n',s1)
                            stop
                    end
                end

                fclose(fid); 
                
                %SAVE the last parameters
                
                LME_param(M,1)=gamma_lme;
                LME_param(M,2)=target_zero;
                LME_param(M,3)=max(2*eps,TolLag);
                LME_param(M,4)=Nelder;
                LME_param(M,5)=tol_search;
                LME_param(M,6)=prop;
                LME_param(M,7)=gamma_top;
                LME_param(M,8)=nb_grade;

        end

        function MAT_POINT=separation(MAT_POINT,ph,SEP,NODE_LIST)
            
            global GEOMETRY
            
            rbs=NODE_LIST.rbs;
            RB=NODE_LIST.RB;
            
            l=SEP(ph).l;
            h=SEP(ph).h*6;
            
            if l~=0
                
                elements=GEOMETRY.mat_points;
                nodes=GEOMETRY.nodes;
                x_a=GEOMETRY.x_0;
                
                for s=1:length(l)
    
                    if l<=rbs
                        rlist=RB{l(s)};
                        x0=[rlist(1),rlist(2)];
                        x1=[rlist(3),rlist(4)];
                        M=(x0(2)-x1(2))/(x0(1)-x1(1));
                        if isinf(M)
                            M1=0;
                        else
                            M1=-1/M;
                        end
                        
                        rM=sqrt(M^2+1);
                        rM1=sqrt(M1^2+1);
                        
                        %Allocate
                        d=zeros(nodes,1);
                        d1=zeros(nodes,1);
                        d2=zeros(nodes,1);
                        clas=zeros(nodes,1);
                        dg=zeros(elements,1);
                        dg1=zeros(elements,1);
                        dg2=zeros(elements,1);
                        clasg=zeros(elements,1);
                        
                        for i=1:nodes
                            % Identify side of the nodes
                            if isinf(M) || abs(M)>1e2 
                                d(i)=x_a(i,1)-x0(1);
                            else
                                
                                d(i)=(M*x_a(i,1)-x_a(i,2)+x0(2)-M*x0(1))/rM;
                            end
                            
                            if isinf(M1) || abs(M1)>1e2 
                                d1(i)=x0(1)-x_a(i,1);
                                d2(i)=x1(1)-x_a(i,1);
                            else
                                d1(i)=(M1*x0(1)-x0(2)+x_a(i,2)-M1*x_a(i,1))/rM1;
                                d2(i)=(M1*x1(1)-x1(2)+x_a(i,2)-M1*x_a(i,1))/rM1;
                            end
                            
                            if d1(i)==0 || d2(i)==0
                                in=1;
                            elseif sign(d1(i))~= sign(d2(i))
                                in=1;
                            else
                                in=0;
                            end
                            
                            if abs(d(i))>h/2 || in==0
                                clas(i)=0;
                            elseif d(i)>=0
                                clas(i)=1;
                            else
                                clas(i)=2;
                            end
                        end
                        
                        for i=1:elements
                            
                            x_g=MAT_POINT{ph}(i).xg;
        
                            % Identify side of the material point
                            if isinf(M) || abs(M)>1e2 
                                dg(i)=x_g(1)-x0(1);
                            else
                                dg(i)=(M*x_g(1)-x_g(2)+x0(2)-M*x0(1))/rM;
                            end
                            
                            if isinf(M1) || abs(M1)>1e2 
                                dg1(i)=x0(1)-x_g(1);
                                dg2(i)=x1(1)-x_g(1);
                            else
                                dg1(i)=(M1*x0(1)-x0(2)+x_g(2)-M1*x_g(1))/rM1;
                                dg2(i)=(M1*x1(1)-x1(2)+x_g(2)-M1*x_g(1))/rM1;
                            end
                            
                            if dg1(i)==0 || dg2(i)==0
                                in=1;
                            elseif sign(dg1(i))~= sign(dg2(i))
                                in=1;
                            else
                                in=0;
                            end
                            
                            if abs(dg(i))>h/2 || in==0
                                clasg(i)=0;
                            elseif dg(i)>=0
                                clasg(i)=1;
                            else
                                clasg(i)=2;
                            end

                            % Initialize list of incompatible nodes
                            no=MAT_POINT{ph}(i).near_p;
                            t=length(no);
                            if no==0
                                t=0;
                            end
                            % Exclude nodes in the other side
                            if clasg(i)~=0
                                for j=1:nodes
                                    if clas(j)~=0 && clas(j)~=clasg(i)
                                        l=0;
                                        for k=1:t
                                            if j==no(k)
                                                l=1;
                                            end
                                        end
                                        if l==0
                                            t=t+1;
                                            no(t)=j;
                                        end
                                    end
                                end
                                MAT_POINT{ph}(i).near_p=no;
                                %plot_nb(i,no,x_a,GEOMETRY.elem_c,1,0);
                                %i;
                            end

                            clear no
                        end

                    else
                        error('Unrecognized Rigid Body')
                    end
                end
            end
            
        end
        
        function MAT_POINT=initialize(MAT_POINT,NODE_LIST)
            
            global GEOMETRY LME_param SOLVER
            
            SEP=LME.read;
            
            [phases,~]=size(SOLVER.PHASES);
            [shf,~]=size(LME_param);
            for sh=1:shf
                for ph=1:phases
                    if SOLVER.PHASES{ph,2}==sh
                        break;
                    end
                end
                
                MAT_POINT=LME.separation(MAT_POINT,ph,SEP,NODE_LIST);
                
                gamma_lme = LME_param(sh,1);
                gamma_=gamma_lme*ones(GEOMETRY.mat_points,1);

                n_sp=LME.nodalspacing(MAT_POINT{ph});

                lam_LME=zeros(GEOMETRY.mat_points,GEOMETRY.sp);
                for i=1:GEOMETRY.mat_points
                    beta_=gamma_(i)/n_sp(i)^2;
                    [lam_LME(i,:)]=LME.first_lambda(i,MAT_POINT{ph}(i).near,...
                        GEOMETRY.xg_0,GEOMETRY.x_0,beta_);%First lambda
                    MAT_POINT{ph}(i).w=gamma_(i);
                    MAT_POINT{ph}(i).xi=lam_LME(i,:);
                end  
                
                for ph1=ph+1:phases
                    if SOLVER.PHASES{ph1,2}==sh
                        for i=1:GEOMETRY.mat_points
                            MAT_POINT{ph1}(i).w = MAT_POINT{ph}(i).w;
                            MAT_POINT{ph1}(i).xi= MAT_POINT{ph}(i).xi;
                        end
                    end
                end
            end
        end
         
        function [MAT_POINT,REMAP]=REMAPPING(MAT_POINT,i)

            global LME_param

            REMAP=0;

            %Tolerance to make new re-mapping or not and proportion to do it
            tol_search = LME_param(5);
            prop = LME_param(6);
            %Minimum gamma
            gamma_top=LME_param(7);

            dis(3)=0;
            gamma_=MAT_POINT(i).w;
            Ep=MAT_POINT(i).EP;
            num=min(Ep(:,1),Ep(:,2));
            for j=1:3
                dis(j)=(Ep(j,1)-Ep(j,2))/num(j);
            end
            Dis=max(abs(dis));
            if Dis>tol_search
                m_dis=max(dis);
                for j=1:3
                    MAT_POINT(i).EP(j,2)=MAT_POINT(i).EP(j,1);
                end
                REMAP=1;
                if m_dis>0
                    gamma_=max(gamma_-m_dis/prop,gamma_top);
                end
            end 
            MAT_POINT(i).w=gamma_;

        end
        
        function h=nodalspacing(MAT_POINT)
            global GEOMETRY
            
            h=zeros(GEOMETRY.mat_points,1);
            
            for i=1:GEOMETRY.mat_points
                x0=MAT_POINT(i).xg;
                e=MAT_POINT(i).element;
                near=[e, GEOMETRY.element_near{e}];
                d=zeros(length(near),1);
                k=0;
                for j=1:length(near)
                    for n=1:GEOMETRY.mat_points
                        if near(j)==MAT_POINT(n).element && n~=i
                            k=k+1;
                            xj=MAT_POINT(n).xg;
                            d(k)=sqrt((x0-xj)*(x0'-xj'));
                        end
                    end
                end
                h(i)=min(d);
                clear d
            end
            
        end
    end
 end


