
function [Shape_function]=LME_EP_ini(Mat_state,Disp_field)
    
    global GEOMETRY SOLVER MATERIAL LME_param
    
    jacobians=Mat_state.J;
    xg=Mat_state.xg;
    
    x_a=Disp_field.x_a;
    
    h=GEOMETRY.h_ini.*sqrt(jacobians);
    volume=GEOMETRY.Area.*jacobians;
    
    B={0};
    near={0};
    [near]=first_near(GEOMETRY.elem,near);
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
    
    %% WRAP or not?
    beta_=zeros(GEOMETRY.elements,1);
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
        [p,dp,lam_LME,FAIL]=shapef...
            (x_a,beta_,xg,near,TolLag,Nelder,wrap,p,lam_LME);
        if FAIL
            fprintf('FAIL in the initial calculation of Shape function \n');
            stop;
        end
    end

    %% MAKE Bbar - Patch
    if SOLVER.B_BAR==1
        [B,near,p]=bbar_matrix...
            (x_a,GEOMETRY.patch_con,GEOMETRY.patch_el,...
            GEOMETRY.elements,near,dp,p,volume,B,wrap,SOLVER.B_BAR);
    else
        [B]=b_m(GEOMETRY.elements,near,dp,GEOMETRY.sp,wrap,B,xg,p);
    end
    
    %% SAVE info
    Shape_function.title='Shape function, derivatives and connectivity';
    Shape_function.B=B;
    Shape_function.near=near;
    Shape_function.p=p;
    Shape_function.gamma=gamma_;
    Shape_function.lambda=lam_LME;
    Shape_function.EP=EP;

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
