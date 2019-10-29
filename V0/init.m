function [ste,ste_p,MAT_POINT,Disp_field,Int_var,Mat_state,GLOBAL,OUTPUT,...
        stiff_mtx,load_s,MATRIX]=init(MAT_POINT)
    
    global GEOMETRY VARIABLE MATERIAL TIME SOLVER
        
    l0          = GEOMETRY.df*GEOMETRY.nodes;
    l1          = GEOMETRY.f_dim;
    l2          = GEOMETRY.s_dim;
    elements    = GEOMETRY.mat_points;
    nodes       = GEOMETRY.nodes;
    dim         = SOLVER.dim; 

    %----------------------------------------------------------------------
    % DEFINITION OF STRUCTURES
    %----------------------------------------------------------------------
    
    Disp_field.title='Displacement field, velocity, acceleration and position';
    Disp_field.d  = zeros(l0,2);
    Disp_field.v  = zeros(l0,2);
    Disp_field.a  = zeros(l0,2);
    Disp_field.x_a= zeros(nodes,2);
    
    Int_var.title='Internal variables: plastic information';
    Int_var.dgamma  = zeros(elements,2);
    Int_var.gamma  = zeros(elements,2);
    Int_var.epsv  = zeros(elements,2);
    Int_var.Sy  = zeros(elements,2);
    Int_var.Sy_r= zeros(elements,2);
    Int_var.H   = zeros(elements,2);
    Int_var.eta = zeros(elements,2);
    Int_var.P0  = zeros(elements,1);
       
    Mat_state.title='Material state: material point information';

    Mat_state.F=zeros(l1*elements,2);
    Mat_state.Be=zeros(l1*elements,2);
    Mat_state.Sigma=zeros(l2*elements,2);
    Mat_state.fint=zeros(l0,2);
    
    if SOLVER.UW==1
        Mat_state.k=zeros(elements,1);
        Mat_state.pw=zeros(elements,2);
        Mat_state.Fw=zeros(l1*elements,2);
    elseif SOLVER.UW==2
        Mat_state.k=zeros(elements,1);
        Mat_state.pw=zeros(elements,2);
        Mat_state.dpw=zeros(GEOMETRY.sp*elements,2);
    end
    
    GLOBAL.title='Variables to save';
    GLOBAL.d  = zeros(l0,dim);
    GLOBAL.v  = zeros(l0,dim);
    GLOBAL.a  = zeros(l0,dim);
    
    GLOBAL.dgamma       = zeros(elements,dim);
    GLOBAL.gamma        = zeros(elements,dim);
    GLOBAL.epsv         = zeros(elements,dim);
    GLOBAL.gamma_nds    = zeros(nodes,dim);
    GLOBAL.Sy           = zeros(elements,dim);
    GLOBAL.Sy_r         = zeros(elements,dim);
    GLOBAL.H            = zeros(elements,dim);
    GLOBAL.eta          = zeros(elements,dim);
    
    GLOBAL.Ps           = zeros(elements,dim);
    GLOBAL.Qs           = zeros(elements,dim);
       
    GLOBAL.F    = zeros(l1*elements,dim);
    GLOBAL.Be   = zeros(l1*elements,dim);
    GLOBAL.Sigma= zeros(l2*elements,dim);
    GLOBAL.Es   = zeros(l2*elements,dim);
    GLOBAL.Es_p = zeros(l2*elements,dim);
    GLOBAL.J    = ones(elements,dim);
    GLOBAL.xg   = zeros(elements*GEOMETRY.sp,dim);
    if SOLVER.UW==1
        GLOBAL.pw=zeros(elements,dim);
        GLOBAL.Fw=zeros(l1*elements,dim);
    end
    
    MATRIX=DYN_MATRIX;
  
    %----------------------------------------------------------------------
    % VECTORS OF PARAMETERS
    %----------------------------------------------------------------------

    for i=1:elements
        Int_var.Sy(i,2) = MATERIAL.MAT(7,MATERIAL.e(i));
        if SOLVER.UW>0
            Mat_state.k(i) = ...
                MATERIAL.MAT(15,MATERIAL.e(i))/VARIABLE.rho_w/VARIABLE.g;
        end
    end
    
    %----------------------------------------------------------------------
    % VECTORS OF ZEROS and ONES   or taken from FILE
    %----------------------------------------------------------------------
    
    if SOLVER.INIT_file~=0
        
        time_aux=TIME.t;     
        
        s=strcat('load("',SOLVER.INIT_file,...
            '.mat","GLOBAL","TIME","ste","ste_p","OUTPUT")');
        eval(s);
        

        if SOLVER.INIT_STEP~=0
            l1=length(TIME.t);
            ste_p=SOLVER.INIT_STEP;
            [~,ste]=min(abs(TIME.t(:)-TIME.tp(ste_p)*ones(l1,1)));
        end
        
        [TIME.t,ste]=fix_time(ste,TIME.t,time_aux,ste_p);
        
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
        Int_var.P0(:,1)=GLOBAL.Ps(:,1);
        
        Mat_state.Sigma(:,2)=GLOBAL.Sigma(:,ste_p);                %STRESS
        Mat_state.F(:,2)=GLOBAL.F(:,ste_p);          %DEFORMATION GRADIENT
        Mat_state.F(:,1)=GLOBAL.F(:,ste_p);          %DEFORMATION GRADIENT
        Mat_state.Be(:,2)=GLOBAL.Be(:,ste_p);              %LEFT CAUCHY GREEN 
        
        [MAT_POINT]=LIB.list2S(MAT_POINT,'J',GLOBAL.J(:,ste_p));
        
        [MAT_POINT]=LIB.list2S(MAT_POINT,'xg',...
            reshape(GLOBAL.xg(:,ste_p),[elements,GEOMETRY.sp]));
        
        if SOLVER.UW==1
            Mat_state.pw(:,2)=GLOBAL.pw(:,ste_p);
            Mat_state.Fw(:,2)=GLOBAL.Fw(:,ste_p);
            Mat_state.Fw(:,1)=GLOBAL.Fw(:,ste_p);
        end
        
        [OUTPUT_2]=read_output;     
        OUTPUT=together_outputs(OUTPUT,OUTPUT_2,dim);
        
        
        MAT_POINT=shape_function_calculation(0,MAT_POINT,Disp_field);
        
        [stiff_mtx,Int_var,Mat_state]=...
        Constitutive(1,ste,Int_var,Mat_state,MAT_POINT);
        [load_s(:,1),~]=calculate_forces...
            (ste,MAT_POINT,Disp_field,Mat_state,OUTPUT,MATRIX);
        SOLVER.INIT_file=0;
        
        Mat_state.fint(:,2)=Mat_state.fint(:,1);
        if SOLVER.UW==1
            Mat_state.pw(:,3)=Mat_state.pw(:,2)-Mat_state.pw(:,1);            
        end
        
    else
        Disp_field.x_a=GEOMETRY.x_0;
        GLOBAL.xg(:,1)=reshape(GEOMETRY.xg_0,[GEOMETRY.sp*GEOMETRY.mat_points,1]);
        
        [Mat_state]=ini_F(Mat_state,l1,elements,GEOMETRY.sp,SOLVER.UW);
        
        TIME.tp(dim,1)=0;
        
        ste=1;
        ste_p=1;
        TIME.tp(1)=TIME.t(ste);
        
        % OUTPUT
        %--------------------------------------------------------------------------
        [OUTPUT]=read_output;
        
        [~,GLOBAL.d(:,1),GLOBAL.v(:,1),~]=calculate_boundaries(1,0);
        
        Disp_field.v(:,2)=GLOBAL.v(:,1);
        Disp_field.a(:,2)=GLOBAL.a(:,1);
        Disp_field.d(:,2)=GLOBAL.d(:,1);
        
        [stiff_mtx,GLOBAL,Disp_field,Int_var,Mat_state,load_s,OUTPUT]=...
            initial_state(MAT_POINT,Disp_field,Int_var,Mat_state,...
            GLOBAL,OUTPUT,MATRIX);
        
    end 

end

function [Mat_state]=ini_F(Mat_state,l1,elements,sp,UW)

    for e=1:elements       
        Mat_state.F((e-1)*l1+1,1)=1; 
        Mat_state.F((e-1)*l1+sp+2,1)=1;
        Mat_state.F(e*l1,1)=1;
        
        Mat_state.Be((e-1)*l1+1,1)=1; 
        Mat_state.Be((e-1)*l1+sp+2,1)=1;
        Mat_state.Be(e*l1,1)=1;  
        
    end
    Mat_state.F(:,2)=Mat_state.F(:,1);
    Mat_state.Be(:,2)=Mat_state.Be(:,1);
    
    if UW==1
        for e=1:elements       
            Mat_state.Fw((e-1)*l1+1,1)=1; 
            Mat_state.Fw((e-1)*l1+sp+2,1)=1;
            Mat_state.Fw(e*l1,1,1)=1;
        end
        Mat_state.Fw(:,2)=Mat_state.Fw(:,1);
    end
end

function [time_2,ste]=fix_time(ste,time,time_aux,ste_p)

    global SOLVER BOUNDARY LOAD
    time_now = time(ste);
    
    l1=length(time_aux);
    [~,j]=min(abs(time_aux(:)-time_now*ones(l1,1)));
    
    if time_aux(j)<=time_now
        ste_2=j+1;
    else
        ste_2=j;
    end
    l2=ste+l1-ste_2+1;
    time_2=zeros(l2,1);
    time_2(1:ste)=time(1:ste);
    time_2(ste+1:l2)=time_aux(ste_2:l1);

    SOLVER.step_final=l2;
    
    [i,b_s]=size(BOUNDARY.b_mult);
    [~,l_s]=size(LOAD.load_mult);
    rest = i - SOLVER.step_final;
    
    if rest>=0
        BOUNDARY.b_mult(1:rest,:) = [];
        LOAD.load_mult(1:rest,:) = [];
    else
        b_mult_2=string(zeros(SOLVER.step_final,b_s));
        load_mult_2=zeros(SOLVER.step_final,l_s);
        b_mult_2(end-i+1:end,:)=BOUNDARY.b_mult(:,:);
        load_mult_2(end-i+1:end,:)=LOAD.load_mult(:,:);
        BOUNDARY.b_mult=b_mult_2;
        LOAD.load_mult=load_mult_2;
    end
    
    dim_2=floor((l1-ste_2)/SOLVER.SAVE_I);
    SOLVER.dim=dim_2+ste_p;
    fprintf('%i plot steps\n',SOLVER.dim);
    fprintf('Save %i times before the final\n',round(SOLVER.dim/SOLVER.SAVE_F));
end

function [OUTPUT3]=together_outputs(OUTPUT,OUTPUT2,dim)

    OUTPUT3=OUTPUT;
    OUTPUT3.name=OUTPUT2.name;

    new_dim=OUTPUT2.number;
    old_dim=OUTPUT.number;
    
    for i=1:new_dim
        if OUTPUT2.type(i,1)==1
            t=0;
            for j=1:old_dim
                if OUTPUT.type(j,1)==1 && OUTPUT2.type(i,2)==OUTPUT.type(j,2)
                    t=1;
                    break;
                end
            end
            if t==0
                OUTPUT3.number=OUTPUT3.number+1;
                OUTPUT3.type(OUTPUT3.number,:)=OUTPUT2.type(i,:);
            end
        elseif OUTPUT2.type(i,1)==2
            t=0;
            for j=1:old_dim
                if OUTPUT.type(j,1)==2
                    list=abs(OUTPUT.ref_list-OUTPUT2.ref_list);
                    if norm(list)<1e-5
                        t=1;
                        break;
                    end
                end
            end
            if t==0
                OUTPUT3.number=OUTPUT3.number+1;
                OUTPUT3.type(OUTPUT3.number,:)=OUTPUT2.type(i,:);
                OUTPUT3.ref_list(OUTPUT3.number,:)=OUTPUT2.ref_list(i,:);
            end
        end        
    end
    
    if OUTPUT3.number~=OUTPUT.number
        OUTPUT3.list(dim,OUTPUT3.number)=0;
    end

end


function [stiff_mtx,GLOBAL,Disp_field,Int_var,Mat_state,load_s,OUTPUT]=...
         initial_state(MAT_POINT,Disp_field,Int_var,Mat_state,...
         GLOBAL,OUTPUT,MATRIX)

    global GEOMETRY SOLVER
    
    % A. Stress and deformation gradient
    [Mat_state,matrix,Int_var,MAT_POINT]=initial_constitutive...
        (MAT_POINT,Mat_state,Int_var);
    
    % B. External forces
    load_s=zeros(GEOMETRY.nodes*GEOMETRY.df,2);
    [load_s(:,1),OUTPUT]=...
        calculate_forces(1,MAT_POINT,Disp_field,Mat_state,OUTPUT,MATRIX);
    
    % C. Initial calculation if required
    if SOLVER.step0==1 || any(load_s(:,1))
        % 1. Calculate forces
        [Mat_state]=internal_forces(MAT_POINT,Mat_state);

        %GT= load_s(:,1)-Mat_state.fint(:,1);
        %[InvK,GT]=apply_conditions(3,1,matrix,GT);

        % 2. Newton-Raphson
        [Disp_field,Mat_state,MAT_POINT]=Newton_Raphson_solver...
            (1,matrix,0,0,load_s,MAT_POINT,Disp_field,Int_var,Mat_state);

        % 3. Calculate initial update of F   
        %[Mat_state,MAT_POINT]=...
            %update_F(Disp_field.d,Mat_state,MAT_POINT);
    end
    
    SOLVER.step0=0;
    % D. Calculate initial stiffness matrix    
    [Mat_state,MAT_POINT]=...
        update_F(Disp_field.d,Mat_state,MAT_POINT);

    [stiff_mtx,Int_var,Mat_state]=...
                Constitutive(2,1,Int_var,Mat_state,MAT_POINT);
     
     % E. Save information
     for e=1:GEOMETRY.mat_points
         [GLOBAL.Ps(e,1),GLOBAL.Qs(e,1)]=LIB.invar(Mat_state.Sigma(:,1),e);   
     end         %PRESSURE
     
     Mat_state.fint(:,2)= Mat_state.fint(:,1);
     Mat_state.Be(:,2)  = Mat_state.Be(:,1);
     GLOBAL.Be(:,1)     = Mat_state.Be(:,1);
     GLOBAL.Sigma(:,1)  = Mat_state.Sigma(:,1);                %STRESS
     Int_var.Sy(:,2)    = Int_var.Sy(:,1);
     GLOBAL.Sy(:,1)     = Int_var.Sy(:,1);
     Int_var.Sy_r(:,2)  = Int_var.Sy_r(:,1);
     GLOBAL.Sy_r(:,1)   = Int_var.Sy_r(:,1);
     Mat_state.F(:,2)   = Mat_state.F(:,1);
     GLOBAL.F(:,1)      = Mat_state.F(:,1);
     
     if SOLVER.UW==1
         GLOBAL.pw(:,1)=Mat_state.pw(:,1);
     end
     
     % F. Plot initial state
     plot_initial_state(GLOBAL.Ps(:,1),Disp_field.d(:,1));
     
     % G. Displacement?
     if SOLVER.INITIAL_d==0
        Disp_field.d=zeros(GEOMETRY.nodes*GEOMETRY.df,2);
     else
        Disp_field.d(:,2)   = Disp_field.d(:,1);
        Disp_field.d(:,1)   = zeros(GEOMETRY.nodes*GEOMETRY.df,1);
        GLOBAL.d(:,1)       = Disp_field.d(:,2);
     end

end

function plot_initial_state(P0,d0)

    global GEOMETRY
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    x_0=GEOMETRY.x_0;
    
    if df==4
        for j=1:GEOMETRY.nodes
            d((j-1)*sp+1:j*sp,1)=d0((j-1)*df+1:(j-1)*df+sp,1);
        end
    else
        d=d0;
    end
    
    Dx=(max(GEOMETRY.x_0(:,1))-min(GEOMETRY.x_0(:,1)))/100;
    Dy=(max(GEOMETRY.x_0(:,2))-min(GEOMETRY.x_0(:,2)))/100;

    [x1,y1]=meshgrid(min(x_0(:,1)):Dx:max(x_0(:,1)),0:Dy:max(x_0(:,2)));
    figure(1)
    
    subplot(121)
    PW1=griddata(GEOMETRY.xg_0(:,1),GEOMETRY.xg_0(:,2),P0,x1,y1,'cubic');
    surf(x1,y1,PW1)
    colorbar
    axis([min(x_0(:,1)),max(x_0(:,1)),min(x_0(:,2)),max(x_0(:,2))])
        

    % Axis dimension  
    x_1=GEOMETRY.x_0;
    for i=1:GEOMETRY.nodes
        x_1(i,1)=x_1(i,1)+d(i*GEOMETRY.sp-1,1);
        x_1(i,2)=x_1(i,2)+d(i*GEOMETRY.sp,1);
    end
    xmax = max(max(GEOMETRY.x_0(:,1)), max(x_1(:,1)));
    ymax = max(max(GEOMETRY.x_0(:,2)), max(x_1(:,2)));
    elaxis = [0 xmax 0 ymax];
    
    % Structure
    subplot(122)
    
    plot_nb(0,0,GEOMETRY.x_0,GEOMETRY.xg_0,GEOMETRY.elem,d,0)
    
    hold on
    plot_nb(0,0,GEOMETRY.x_0,GEOMETRY.xg_0,GEOMETRY.elem,d,1)
    axis(elaxis)
    title('Current and deformed configuration')
    hold off
    drawnow

end



