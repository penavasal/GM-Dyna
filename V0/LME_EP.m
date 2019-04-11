
function [Shape_function]=LME_EP(Mat_state,Disp_field,INITIAL)
    
    global GEOMETRY SOLVER MATERIAL LME_param
    
    jacobians=Mat_state.J;
    xg=Mat_state.xg;
    
    x_a=Disp_field.x_a;
    
    h=GEOMETRY.h_ini.*sqrt(jacobians);
    volume=GEOMETRY.Area.*jacobians;
    
    B={0};
    near={0};
    p={0};

    % As described in [1], gamma measures the degree of locality, the larger,
    % the closer to Delaunay finite elements
    gamma_lme = LME_param(1);
    gamma_=gamma_lme*ones(GEOMETRY.elements,1);
        
    % This sets the numerical threshold for the support of the shape functions
    target_zero= LME_param(2); 
    % This is the tolerance for the Newton iteration to compute the shape
    % functions
    TolLag= LME_param(3);
    % Nelder Mead method to obtain minimum lambda in shape functions
    Nelder=LME_param(4);  
    
    %Tolerance to make new re-mapping or not and proportion to do it
    tol_search = LME_param(5);
    prop = LME_param(6);
    %Minimum gamma
    gamma_top=LME_param(7);
    
    %% WRAP or not?
    beta_=zeros(GEOMETRY.elements,1);
    if INITIAL==0
        wrap=zeros(GEOMETRY.elements,1);
        EP=zeros(3*GEOMETRY.elements,2);
        for i=1:GEOMETRY.elements
            wrap(i)=1;
            for j=1:3
                EP((i-1)*3+j,1)=1;
                EP((i-1)*3+j,2)=1;
            end  
            beta_(i)=gamma_(i)/h(i)^2;
        end   
        [lam_LME]=first_lambda(GEOMETRY.elem,xg,x_a,beta_);   %First lambda
    else
        dis(3)=0;
        for i=1:GEOMETRY.elements
            if wrap(i)==1
                REMAP=1;
            else
                Ep=EP((i-1)*3+1:(i-1)*3+3,:);
                num=min(Ep(:,1),Ep(:,2));
                for j=1:3
                    dis(j)=(Ep(j,1)-Ep(j,2))/num(j);
                end
                Dis=max(abs(dis));
                if Dis>tol_search
                    m_dis=max(dis);
                    for j=1:3
                        EP((i-1)*3+j,2)=EP((i-1)*3+j,1);
                    end
                    wrap(i)=1;
                    %plot_nb(i, near, x_a, xg, elem_0)
                    REMAP=1;
                    if m_dis>0
                        gamma_(i)=max(gamma_(i)-m_dis/prop,gamma_top);
                    end
                end 
            end
            beta_(i)=gamma_(i)/h(i)^2;
        end
    end
    
    %% MAKE NEAR
    range1=zeros(GEOMETRY.elements,1);
    range_lme=zeros(GEOMETRY.elements,1);
    for i=1:GEOMETRY.elements
        range1(i)=h(i)*sqrt(-1/gamma_(i)*log(target_zero));
        range_lme(i)=max(range1(i)*1, 2*h(i));
    end
    [near]=make_near(x_a,xg,range_lme,wrap,near,MATERIAL.e,MATERIAL.n);


    %% MAKE SHAPE FUNCTION
    if GEOMETRY.sp==2
        [p,dp,lam_LME,FAIL]=LME_shapef...
            (x_a,beta_,xg,near,TolLag,Nelder,wrap,p,lam_LME);
        if FAIL
            fprintf('FAIL in the initial calculation of Shape function \n');
            stop;
        end
    end

end

function [near]=make_near(x_a,x_sample,range,wrap,near,Mat,Mat_nds)
    % Naive construction of the neighbor list
    [n_a,ndim]=size(x_a);
    [n_sample,~]=size(x_sample);

    for i=1:n_sample
        if wrap(i)==1
              near_=near{i};
              mm=length(near_);
              clear near(i)
              if (ndim==1)
                x=x_sample(i);
              else
                x=x_sample(i,:);
              end
              %Naive find in range
              dist=zeros(n_a,1);
              for id=1:ndim
                dist=dist + (x_a(:,id)-x(id)).^2;
              end
              dist=sqrt(dist);
              
              for j=1:n_a
                  %if dist(j)<range(i)
                  if (Mat(i)==Mat_nds(j)) && (dist(j)<range(i))
                      t=0;
                      for k=1:mm
                          if near_(k)==j
                              t=1;
                          end
                      end
                      if t==0;
                          mm=mm+1;
                          near_(mm)=j;
                      end
                  end
              end
            near(i)={near_};
        end
    end
    clear dist
end

function [p,dp,lambda,FAIL]=LME_shapef...
    (x_a,beta_,x_sample,near_,TolLag,Nelder,wrap,p,lambda)
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


    [n_sample,ndim]=size(x_sample);
    FAIL=0;

    for i=1:n_sample        
        if wrap(i)==1           
            clear p(i) p_a
        
              if (ndim==1)
                x=x_sample(i);
              else
                x=x_sample(i,:);
              end
              beta=beta_(i);
              near=near_{i};

              lama=lambda(i*ndim-1:i*ndim);
              lam=lama';

              if Nelder
                  [p_a,J,cnt,lam,tol]=Nelder_LME_mod(x,x_a,beta,near,TolLag,lam,ndim);
              else
                  [p_a,J,~]=NR_LME(x,x_a,beta,near,TolLag,lam,ndim);
              end
              
              if cnt>=250
                  fprintf('Convergence problem in element %i with tol %i \n'...
                      ,i,tol);
              end
              
              lambda(i*ndim-1:i*ndim)=lam';
              p(i)={p_a}; 

              if isnan(p_a)
                  fprintf('Not a number in element %i \n',i);
                  FAIL=1;
                  break;
              elseif isnan(J)
                  fprintf('Not a number in element %i \n',i);
                  FAIL=1;
                  break;
              else
                    %Spacial Gradients
                  dp(i)={ zeros(length(near),ndim) };
                  for ia=1:length(near)
                    if ndim==1
                      dp{i}(ia,1)= -(p_a(ia) * (x-x_a(near(ia))') )/J(1,1);      
                    else
                      dp{i}(ia,:)= -J\(p_a(ia) * (x(:)-x_a(near(ia),:)') );
                    end
                  end
              end
        else
        	dp(i)={0};
        end
    end   %nsample points
    
end

function [near]=first_near(elem,near)
    [n_sample,~]=size(elem);
    for i=1:n_sample
        near(i)={elem(i,:)};
    end
end

function [lambda]=first_lambda(elem,x_sample,x_a,beta_)

[n_sample,bg]=size(elem);
[~,sp]=size(x_a);

lambda=zeros(sp*n_sample,1);

    for i=1:n_sample
        
          x=x_sample(i,:);
          beta=beta_(i);

        % Initialize lambda
          y=zeros(bg,sp);
          for j=1:bg
            for k=1:sp
                y(j,k)=x(k)-x_a(elem(i,j),k);
            end
          end
          A(1,1)=y(2,1)-y(1,1);
          A(1,2)=y(2,2)-y(1,2);
          A(2,1)=y(3,1)-y(1,1);
          A(2,2)=y(3,2)-y(1,2);
          B(1,1)=-beta*(y(1,1)^2+y(1,2)^2-(y(2,1)^2+y(2,2)^2));
          B(2,1)=-beta*(y(1,1)^2+y(1,2)^2-(y(3,1)^2+y(3,2)^2));
          lam=A\B;
          
          lambda(i*sp-1:i*sp)=lam;
    end
end

function [p_a,J,niter]=NR_LME(x,x_a,beta,near,TolLag,lam,ndim)


  R=10*ones(1,ndim);
  niter=0;
  I=eye(ndim);
  
  %Newton iteration
  iflag=0;
  dlam=10*ones(1,ndim);
  while (norm(R)>TolLag)
    [gam,R,J,p_a]=Gamma_(ndim,x_a,x,beta,lam,near);
    if (abs(rcond(J)))<1e-8
      iflag=1;
      disp('Newton Failed, near to singular matrix')
    end
    inv=-J\I;
    dlam=-J\R';
    lam=lam+dlam';
    niter=niter+1;
    if (niter>100) 
      i, niter
      iflag=1;
      disp('Newton Failed 2, no convergence in 100 iterations')
    end
  end
  
end
  
function [p_a,J,cnt,lam,TOL]=Nelder_LME_mod(x,x_a,beta,near,TolLag,lam,ndim)
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
        [f(i),~,~,~]=Gamma_(ndim,x_a,x,beta,lam,near);
    end


    % 1- Ordenar

    [LL,f_m,f_mm,f_w]=order(L,f);
    
    %%%%% ROTATION ?? %%%%%
    ROT=1;
    iter=1;
    while ROT==1 && iter<100
        [Vn,diam]=von(LL);
        if Vn < epsilon
            ROT=1;
            [LL,f_m,f_mm,f_w]=rota2D(LL,near,x_a,x,beta,diam);
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

        [f_r,~,~,~]=Gamma_(ndim,x_a,x,beta,x_r,near);

        if f_r<f_m   
            %3-Expansion
            x_e=x_bar+xi*(x_r-x_bar);
            %x_e=x_r+xi*(x_r-x_bar);
            [f_e,~,~,~]=Gamma_(ndim,x_a,x,beta,x_e,near);

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
                [f_c,~,~,~]=Gamma_(ndim,x_a,x,beta,x_c,near);
                if f_c<f_r
                    cont=1;
                end

            else
                %4- Contraccion dentro
                x_c=x_bar-gam*(x_bar-LL(3,:));
                [f_c,~,~,~]=Gamma_(ndim,x_a,x,beta,x_c,near);
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
                    [f(i),~,~,~]=Gamma_(ndim,x_a,x,beta,lam,near);
                end
                f(1)=f_m;
                [LL,f_m,f_mm,f_w]=order(LL,f);
            end
        end
        cnt=cnt+1;
        %%%%% ROTATION ?? %%%%%
        ROT=1;
        iter=1;
        while ROT==1 && iter<100
            [Vn,diam]=von(LL);
            if Vn < epsilon
                ROT=1;
                [LL,f_m,f_mm,f_w]=rota2D(LL,near,x_a,x,beta,diam);
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
    [~,~,J,p_a]=Gamma_(ndim,x_a,x,beta,lam,near);
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
    [f_1(2),~,~,~]=Gamma_(sp,x_a,x,beta,LL1(2,:),near);
    [f_1(3),~,~,~]=Gamma_(sp,x_a,x,beta,LL1(3,:),near);
    [f_2(2),~,~,~]=Gamma_(sp,x_a,x,beta,LL2(2,:),near);
    [f_2(3),~,~,~]=Gamma_(sp,x_a,x,beta,LL2(3,:),near);
    
    f1=min(f_1);
    f2=min(f_2);
    
    if f1<f2
        LL2=LL1;
        f_2=f_1;
    end
    [f_2(1),~,~,~]=Gamma_(sp,x_a,x,beta,LL1(1,:),near);
    
    [LL2,f_m,f_mm,f_w]=order(LL2,f_2);
end



