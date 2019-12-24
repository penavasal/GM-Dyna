
function [ste,ste_p,MAT_POINT,GLOBAL]=Explicit_solver(...
            BLCK,ste,ste_p,MAT_POINT,Disp_field,Int_var,...
                Mat_state,GLOBAL,MATRIX,load_s)

    %--------------------------------------------------------------------------
    % Initialize variables
    %--------------------------------------------------------------------------
    global SOLVER TIME
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ste=ste+1:SOLVER.step_final(BLCK)
        
        if ste==1
            time_step=TIME{BLCK}.t(ste+1)-TIME{BLCK}.t(ste);
        else
            time_step=TIME{BLCK}.t(ste)-TIME{BLCK}.t(ste-1);
        end 
        gamma=TIME{BLCK}.gamma;
        

        % 1. Forces
        load_s(:,2)=load_s(:,1);
        [load_s(:,1),out_list1]=calculate_forces...
            (ste,MAT_POINT,Disp_field,Mat_state,MATRIX,BLCK);

        % --------------------------------------------------------
        % 2. Predictor         
        [Disp_field,Mat_state,MAT_POINT]=explicit_predictor...
            (ste,Disp_field,MAT_POINT,Mat_state,time_step,gamma);
        

        % 3. Recompute mass and damping matrices
        [MATRIX] = MATRIX.lumped_mass(MAT_POINT,MATRIX,BLCK);
        [MATRIX] = MATRIX.lumped_damp(MAT_POINT,Mat_state,MATRIX);

        % 4. Constitutive &/O Stiffness_mat
        [~,Int_var,Mat_state]=...
            Constitutive(3,ste,Int_var,Mat_state,MAT_POINT,BLCK);

        % 5. Final conditions: corrector
        
        [Disp_field]=explicit_corrector...
            (ste,MATRIX,Disp_field,load_s,Mat_state.fint,time_step,gamma);

        % 6. Storage
        if rem(ste,SOLVER.SAVE_I)==0 || (ste==SOLVER.step_final(BLCK)) 
            ste_p=ste_p+1;
            fprintf('ste_p %i \n',ste_p);
          
            GLOBAL=VECTORS_Store(GLOBAL,ste_p,BLCK,Mat_state,Disp_field,...
                Int_var,MAT_POINT,out_list1);
        end

        % 7. Save info
        if ((rem(ste/SOLVER.SAVE_I,SOLVER.SAVE_F)==0) ...
                || (ste==SOLVER.step_final) || (SOLVER.FAIL==1))
            save(SOLVER.Output(BLCK),'MAT_POINT','GLOBAL','-append')
        end
        
        % 8. Update
        [Disp_field,Mat_state,Int_var]=VECTORS.Update(...
                Disp_field,Mat_state,Int_var);

    end

end

