
function [STEP,MAT_POINT,GLOBAL]=Explicit_solver(...
            STEP,MAT_POINT,Disp_field,Int_var,Mat_state,GLOBAL,MATRIX,load_s)

    %--------------------------------------------------------------------------
    % Initialize variables
    %--------------------------------------------------------------------------
    global SOLVER
    
    BLCK=STEP.BLCK;
    STEP.dt = STEP.dt*SOLVER.time_factor(BLCK);
    STEP.t  = STEP.t + STEP.dt;
    
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
        [Disp_field,Mat_state,MAT_POINT]=explicit_predictor...
            (STEP,Disp_field,MAT_POINT,Mat_state);
        

        % 3. Recompute mass and damping matrices
        [MATRIX] = MATRIX.lumped_mass(MAT_POINT,MATRIX,BLCK);
        [MATRIX] = MATRIX.lumped_damp(MAT_POINT,Mat_state,MATRIX);

        % 4. Constitutive &/O Stiffness_mat
        [~,Int_var,Mat_state]=...
            Constitutive(3,STEP,Int_var,Mat_state,MAT_POINT);

        % 5. Final conditions: corrector
        
        [Disp_field]=explicit_corrector...
            (STEP,MATRIX,Disp_field,load_s,Mat_state.fint);

        % 6. Storage
        if rem(STEP.ste,SOLVER.SAVE_I)==0
            
            STEP.ste_p=STEP.ste_p+1;
            fprintf('ste_p %i \n',STEP.ste_p);
            
            GLOBAL=VECTORS.Store(GLOBAL,STEP,Mat_state,Disp_field,...
                Int_var,MAT_POINT,out_list1);
        end

        % 7. Save info
        if ((rem(STEP.ste/SOLVER.SAVE_I,SOLVER.SAVE_F)==0) || (SOLVER.FAIL==1))
            save(SOLVER.Output(BLCK),'MAT_POINT','GLOBAL','-append')
            if SOLVER.FAIL==1
                stop
            end  
        end
        
        % 8. Update
        [Disp_field,Mat_state,Int_var]=VECTORS.Update(...
                Disp_field,Mat_state,Int_var);
            
        STEP.dt = STEP.dt*SOLVER.time_factor(BLCK);
        STEP.t  = STEP.t + STEP.dt;
    end
    
    % SAVING LAST VALUES
    STEP.ste_p=STEP.ste_p+1;
    fprintf('ste_p %i \n',STEP.ste_p);

    GLOBAL=VECTORS.Store(GLOBAL,STEP,Mat_state,Disp_field,Int_var,...
        MAT_POINT,out_list1);
    save(SOLVER.Output(BLCK),'MAT_POINT','GLOBAL','-append')

end

