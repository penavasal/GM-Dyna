
function [stiff_mtx,GLOBAL,Disp_field,Int_var,Mat_state,load_s,OUTPUT]=...
         initial_state(Shape_function,Disp_field,Int_var,Mat_state,...
         GLOBAL,OUTPUT)

    global GEOMETRY SOLVER
    
    % 1. Stress and deformation gradient
    [Mat_state,matrix,Int_var]=initial_constitutive...
        (Shape_function,Mat_state,Int_var);
    
    % 2. Calculate forces
    load_s=zeros(GEOMETRY.nodes*GEOMETRY.df,2);
    
    [Mat_state]=internal_forces(Shape_function,Mat_state);
    [load_s(:,1),OUTPUT]=...
            calculate_forces(1,Shape_function,Mat_state,Disp_field,OUTPUT);
    
    GT= load_s(:,1)-Mat_state.fint(:,1);
    [InvK,GT]=apply_conditions(3,1,matrix,GT);
    
    % 3. Newton-Raphson
    [Disp_field,Mat_state,Shape_function,FAIL]=Newton_Raphson_solver...
        (1,GT,InvK,0,0,load_s,Shape_function,Disp_field,Int_var,Mat_state,0);

    
    % 4. Calculate initial stiffness matrix    
    [Mat_state,Shape_function]=...
        update_F(Disp_field.d,Mat_state,Shape_function);

    [stiff_mtx,Int_var,Mat_state,~]=...
                Constitutive(2,1,Int_var,Mat_state,Shape_function,FAIL);
     
     % 5. Save information
     for e=1:GEOMETRY.elements
         [GLOBAL.Ps(e,1),GLOBAL.Qs(e,1)]=invar(Mat_state.Sigma(:,1),e);   
     end         %PRESSURE
     
     Mat_state.Be(:,2)  = Mat_state.Be(:,1);
     GLOBAL.Be(:,1)     = Mat_state.Be(:,1);
     GLOBAL.Sigma(:,1)  = Mat_state.Sigma(:,1);                %STRESS
     Int_var.Sy(:,2)    = Int_var.Sy(:,1);
     GLOBAL.Sy(:,1)     = Int_var.Sy(:,1);
     Int_var.Sy_r(:,2)  = Int_var.Sy_r(:,1);
     GLOBAL.Sy_r(:,1)   = Int_var.Sy_r(:,1);
     Mat_state.F(:,2)   = Mat_state.F(:,1);
     GLOBAL.F(:,1)      = Mat_state.F(:,1);
     
     if SOLVER.UW
         Mat_state.pw(:,3)=SOLVER.INITIAL_COND(2)*...
             ones(GEOMETRY.elements,1)-Mat_state.pw(:,1);
         GLOBAL.pw(:,1)=Mat_state.pw(:,1)+Mat_state.pw(:,3);
     end
     
     % 6. Plot initial state
     plot_initial_state(GLOBAL.Ps(:,1),Disp_field.d(:,1),Shape_function.near);
     
     % 7. Displacement?
     if SOLVER.INITIAL_d==0
        Disp_field.d=zeros(GEOMETRY.nodes*GEOMETRY.df,2);
     else
        Disp_field.d(:,2)   = Disp_field.d(:,1);
        Disp_field.d(:,1)   = zeros(GEOMETRY.nodes*GEOMETRY.df,1);
        GLOBAL.d(:,1)       = Disp_field.d(:,2);
     end

end

function plot_initial_state(P0,d0,near)

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
    PW1=griddata(GEOMETRY.xg0(:,1),GEOMETRY.xg0(:,2),P0,x1,y1,'cubic');
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
    
    plot_nb(0,near,GEOMETRY.x_0,GEOMETRY.xg0,GEOMETRY.elem,d,0)
    
    hold on
    plot_nb(0,near,GEOMETRY.x_0,GEOMETRY.xg0,GEOMETRY.elem,d,1)
    axis(elaxis)
    title('Current and deformed configuration')
    hold off
    drawnow

end