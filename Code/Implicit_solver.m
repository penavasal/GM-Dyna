
function [ste,ste_p,MAT_POINT,GLOBAL]=Implicit_solver(...
            BLCK,ste,ste_p,MAT_POINT,Disp_field,Int_var,...
                Mat_state,GLOBAL,stiff_mtx,MATRIX,load_s)

    %--------------------------------------------------------------------------
    % Initialize variables
    %--------------------------------------------------------------------------
    global SOLVER
    
    %--------------------------------------------------------------------------
    % Initial matrices
    %--------------------------------------------------------------------------
    if SOLVER.DYN(BLCK)==1
        [MATRIX]=MATRIX.matrices(Mat_state,MAT_POINT,Disp_field.d,MATRIX,BLCK);
    end
         
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % LOOP
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ste=ste+1:SOLVER.step_final(BLCK)-1

        % 1. Forces
        load_s(:,2)=load_s(:,1);
        [load_s(:,1),out_list1]=calculate_forces...
            (ste,MAT_POINT,Disp_field,Mat_state,MATRIX,BLCK);
       
        % --------------------------------------------------------
        % 2. Predictor      
        [Disp_field,Mat_state,MAT_POINT]=implicit_predictor...
            (ste,stiff_mtx,MATRIX.mass,MATRIX.damp,load_s,Disp_field,...
            MAT_POINT,Mat_state,Int_var,BLCK);
    
        % --------------------------------------------------------

        % 3. Recompute mass and damping matrices
        if SOLVER.DYN(BLCK)==1
            [MATRIX]=MATRIX.matrices(Mat_state,MAT_POINT,Disp_field.d,MATRIX,BLCK);
        end
        
        % 4. Constitutive & Stiffness_mat
        [stiff_mtx,Int_var,Mat_state]=...
                Constitutive(2,ste,Int_var,Mat_state,MAT_POINT,BLCK);
        
        % 5. Storage
        if rem(ste,SOLVER.SAVE_I)==0 || (ste==SOLVER.step_final(BLCK)-1) 
            ste_p=ste_p+1;
            fprintf('ste_p %i \n',ste_p);
          
            GLOBAL=VECTORS.Store(GLOBAL,ste_p,ste,BLCK,Mat_state,Disp_field,...
                Int_var,MAT_POINT,out_list1);
        end

        % 6. Save info
        if ((rem(ste/SOLVER.SAVE_I,SOLVER.SAVE_F)==0) ...
                || (ste==SOLVER.step_final(BLCK)) || (SOLVER.FAIL==1))
            save(SOLVER.Output(BLCK),'MAT_POINT','GLOBAL','-append')
            if SOLVER.FAIL==1
                stop
            end   
        end
        
        % 7. Update
        [Disp_field,Mat_state,Int_var]=VECTORS.Update(...
                Disp_field,Mat_state,Int_var);

     end
     
end
