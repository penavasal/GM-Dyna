
function Implicit_solver

    tic;
    
    %--------------------------------------------------------------------------
    % Initialize variables
    %--------------------------------------------------------------------------
    global GEOMETRY SOLVER TIME VARIABLE

    %----------------------------------------------------------------------
    % Initial state & initial shape functions and matrixes - Save
    %----------------------------------------------------------------------
    
    [ste,ste_p,Shape_function,Disp_field,Int_var,Mat_state,GLOBAL,OUTPUT,...
        stiff_mtx,load_s]=init;
    
    
    save(OUTPUT.name, 'GEOMETRY', 'VARIABLE', 'SOLVER');

    %--------------------------------------------------------------------------
    % Initial matrices
    %--------------------------------------------------------------------------
    [mass_mtx,damp_mtx]=dyn_matrices(Mat_state,Shape_function,Disp_field.d);
    
    [matrix]=G_matrix(mass_mtx,stiff_mtx,damp_mtx,ste);
        
    [InvK,~]=apply_conditions(0,ste,matrix,0);
         
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % LOOP
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FAIL=0;
    for ste=ste+1:SOLVER.step_final

        % 1. Forces
        load_s(:,2)=load_s(:,1);
        [load_s(:,1),OUTPUT]=calculate_forces...
            (ste,Shape_function,Mat_state,Disp_field,OUTPUT);
        
        [GT]=G_calculation(Disp_field.d,Disp_field.a,Disp_field.v,...
            Mat_state.fint,mass_mtx,damp_mtx,load_s(:,1),load_s(:,2),ste);  

        [~,GT]=apply_conditions(1,ste,matrix,GT);

        % --------------------------------------------------------
        % 2. Predictor      
        [Disp_field,Mat_state,Shape_function,FAIL]=implicit_predictor...
            (ste,GT,InvK,mass_mtx,damp_mtx,load_s,Disp_field,...
            Shape_function,Mat_state,Int_var,FAIL);
    
        % --------------------------------------------------------

        % 3. Recompute mass and damping matrices
        [mass_mtx,damp_mtx]=dyn_matrices(Mat_state,Shape_function,Disp_field.d);

        % 4. Constitutive & Stiffness_mat
        [stiff_mtx,Int_var,Mat_state,~]=...
                Constitutive(2,ste,Int_var,Mat_state,Shape_function,FAIL);
        
        [OUTPUT]=reaction(Mat_state.fint,OUTPUT);
        
        %  5. Assemble time integration matrix and apply conditions
        [matrix]=G_matrix(mass_mtx,stiff_mtx,damp_mtx,ste);
        [InvK,~]=apply_conditions(0,ste,matrix,0);
 
        % 6. Storage
        if rem(ste,SOLVER.SAVE_I)==0

            ste_p=ste_p+1;
            fprintf('ste_p %i \n',ste_p);

            for i=1:OUTPUT.number
                OUTPUT.list(ste_p,i)=OUTPUT.inst(i);
            end
    
            GLOBAL.d(:,ste_p)   = Disp_field.d(:,1);
            GLOBAL.a(:,ste_p)   = Disp_field.a(:,1);
            GLOBAL.v(:,ste_p)   = Disp_field.v(:,1);
            
            GLOBAL.gamma(:,ste_p)   = Int_var.gamma(:,1);
            GLOBAL.dgamma(:,ste_p)  = Int_var.dgamma(:,1);
            GLOBAL.Sy(:,ste_p)      = Int_var.Sy(:,1);
            GLOBAL.Sy_r(:,ste_p)    = Int_var.Sy_r(:,1);

            [GLOBAL.gamma_nds(:,ste_p)]=Ep2Ep_n...
                (GLOBAL.gamma,Shape_function,ste_p);

            for e=1:GEOMETRY.elements
                [GLOBAL.Ps(e,ste_p),GLOBAL.Qs(e,ste_p)]=...
                    invar(Mat_state.Sigma(:,1),e);   %PRESSURE
            end         
            [GLOBAL.Es(:,ste_p),GLOBAL.Es_p(:,ste_p)]=...
                strains(Mat_state.F,Mat_state.Be);
            GLOBAL.J(:,ste_p)       = Mat_state.J(:,1);        %JACOBIAN
            GLOBAL.Sigma(:,ste_p)   = Mat_state.Sigma(:,1);    %STRESS
            GLOBAL.F(:,ste_p)       = Mat_state.F(:,1);  %DEFORMATION GRADIENT
            GLOBAL.Be(:,ste_p)      = Mat_state.Be(:,1); %LEFT CAUCHY GREEN 
            if SOLVER.UW
                GLOBAL.pw(:,ste_p) = Mat_state.pw(:,1)+Mat_state.pw(:,3);
                GLOBAL.Fw(:,ste_p) = Mat_state.Fw(:,1);
            end

            TIME.tp(ste_p,1)=TIME.t(ste);
      
        end

        % 7. Save info
        if ((rem(ste/SOLVER.SAVE_I,SOLVER.SAVE_F)==0) ...
                || (ste==SOLVER.step_final) || (FAIL==1))
            GLOBAL.xg  = Mat_state.xg;
            save(OUTPUT.name,'ste','ste_p','TIME','Shape_function',...
                'GLOBAL','OUTPUT','-append')
        end
        
        % 8.. Update
        Disp_field.d(:,2)=Disp_field.d(:,1);
        Disp_field.v(:,2)=Disp_field.v(:,1);
        Disp_field.a(:,2)=Disp_field.a(:,1);
        
        Mat_state.Sigma(:,2)=Mat_state.Sigma(:,1);                
        Mat_state.F(:,2)=Mat_state.F(:,1);
        Mat_state.Be(:,2)=Mat_state.Be(:,1);
        Mat_state.fint(:,2)=Mat_state.fint(:,1);
        if SOLVER.UW
        	Mat_state.pw(:,2)=Mat_state.pw(:,1);
            Mat_state.Fw(:,2)=Mat_state.Fw(:,1);
        end
        
        Int_var.gamma(:,2)  = Int_var.gamma(:,1);
        Int_var.Sy(:,2)     = Int_var.Sy(:,1);
        Int_var.Sy_r(:,2)   = Int_var.Sy_r(:,1);
        Int_var.dgamma(:,2) = Int_var.dgamma(:,1);

     end
     
     tfin=toc;
     
     disp(tfin);


end

