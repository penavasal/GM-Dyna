
function SOLVER(MAT_POINT)
    
    tic;

    global SOLVER VARIABLE GEOMETRY
     
    %----------------------------------------------------------------------
    % Initial state, matrixes and load
    %----------------------------------------------------------------------    
    [STEP,MAT_POINT,Disp_field,Int_var,Mat_state,GLOBAL,stiff_mtx]=...
        init(MAT_POINT);
    
    %perturbated(MAT_POINT,Mat_state,Disp_field,Int_var,stiff_mtx,STEP);
    
    MATRIX=DYN_MATRIX; 
    
    [load_s,out_list1]=...
        calculate_forces(STEP,MAT_POINT,Mat_state,MATRIX);
    GLOBAL.OutputList(STEP.ste_p,:)=out_list1;
    
    %----------------------------------------------------------------------
    % Solver for each Block
    %---------------------------------------------------------------------- 
    while STEP.BLCK<SOLVER.BLOCKS+1
        
        if STEP.BLCK>1
            if SOLVER.Output(STEP.BLCK)~=SOLVER.Output(STEP.BLCK-1)  || ...
                    isfile(SOLVER.OutputSTEP.(BLCK))==0
                save(SOLVER.Output(STEP.BLCK), 'GEOMETRY', 'VARIABLE', 'SOLVER');
            end
        else
            save(SOLVER.Output(STEP.BLCK), 'GEOMETRY', 'VARIABLE', 'SOLVER');
        end

        if SOLVER.IMPLICIT(STEP.BLCK)==0
            [STEP,MAT_POINT,GLOBAL]=...
                Explicit_solver(STEP,MAT_POINT,Disp_field,Int_var,...
                Mat_state,GLOBAL,MATRIX,load_s);
        else
            [STEP,MAT_POINT,GLOBAL]=...
                Implicit_solver(STEP,MAT_POINT,Disp_field,Int_var,...
                Mat_state,GLOBAL,stiff_mtx,MATRIX,load_s);
        end
        STEP.BLCK=STEP.BLCK+1;
        
        if STEP.BLCK>SOLVER.BLOCKS
            break;
        else
            % Update
            [Disp_field,Mat_state,Int_var,stiff_mtx]=VECTORS.Update_ini(...
                STEP,GLOBAL,Disp_field,Mat_state,Int_var,MAT_POINT);
            
            [load_s,out_list1]=...
                calculate_forces(STEP,MAT_POINT,Mat_state,MATRIX);
            GLOBAL.OutputList(STEP.ste_p,:)=out_list1;
        end
    end
    
    tfin=toc;

    disp(tfin);

end