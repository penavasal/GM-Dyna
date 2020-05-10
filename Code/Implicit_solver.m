
function [STEP,MAT_POINT,GLOBAL]=Implicit_solver(STEP,MAT_POINT,...
            Disp_field,Int_var,Mat_state,GLOBAL,stiff_mtx,MATRIX,load_s)

    %--------------------------------------------------------------------------
    % Initialize variables
    %--------------------------------------------------------------------------
    global SOLVER
    BLCK=STEP.BLCK;
    
    STEP.dt = SOLVER.time_step(BLCK)*SOLVER.time_factor(BLCK);
    STEP.t  = STEP.t + STEP.dt;
    
    %--------------------------------------------------------------------------
    % Initial matrices
    %--------------------------------------------------------------------------
    if SOLVER.DYN(BLCK)==1
        [MATRIX]=MATRIX.matrices(Mat_state,MAT_POINT,Disp_field.d,MATRIX,BLCK);
    end
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while STEP.t < SOLVER.Time_final(BLCK) 
        STEP.ste=STEP.ste+1;

        % 1. Forces
        load_s(:,2)=load_s(:,1);
        [load_s(:,1),out_list1]=calculate_forces...
            (STEP,MAT_POINT,Mat_state,MATRIX);
       
        % --------------------------------------------------------
        % 2. Predictor      
        [Disp_field,Mat_state,MAT_POINT,STEP]=implicit_predictor...
            (STEP,stiff_mtx,MATRIX.mass,MATRIX.damp,load_s,Disp_field,...
            MAT_POINT,Mat_state,Int_var);
    
        % --------------------------------------------------------

        % 3. Recompute mass and damping matrices
        if SOLVER.DYN(BLCK)==1
            [MATRIX]=MATRIX.matrices(Mat_state,MAT_POINT,Disp_field.d,MATRIX,BLCK);
        end
        
        % 4. Constitutive & Stiffness_mat
        [stiff_mtx,Int_var,Mat_state,STEP]=...
                Constitutive.update(2,STEP,Int_var,Mat_state,MAT_POINT);
        
        % 5. Storage
        if rem(STEP.ste,SOLVER.SAVE_I)==0
            STEP.ste_p=STEP.ste_p+1;
            fprintf('ste_p %i \n',STEP.ste_p);
          
            [GLOBAL,STEP]=VECTORS.Store(GLOBAL,STEP,Mat_state,Disp_field,...
                Int_var,MAT_POINT,out_list1);
        end

        % 6. Save info
        if ((rem(STEP.ste/SOLVER.SAVE_I,SOLVER.SAVE_F)==0) || (STEP.FAIL==1))
            save(SOLVER.Output(BLCK),'MAT_POINT','GLOBAL','-append')
            if STEP.FAIL==1
                error('Simulation failed');
            end  
        end
        
        % 7. Update
        STEP=Time_Scheme.step(STEP,Disp_field);
        [Disp_field,Mat_state,Int_var]=VECTORS.Update(...
                Disp_field,Mat_state,Int_var);

    end
     
    % SAVING LAST VALUES
    STEP.ste_p=STEP.ste_p+1;
    fprintf('ste_p %i \n',STEP.ste_p);
    
    GLOBAL.final_block(STEP.BLCK)=GLOBAL.ste_p;

    [GLOBAL,STEP]=VECTORS.Store(GLOBAL,STEP,Mat_state,Disp_field,Int_var,...
        MAT_POINT,out_list1);
    save(SOLVER.Output(BLCK),'MAT_POINT','GLOBAL','-append')
     
end
