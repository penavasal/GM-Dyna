
function SOLVER(MAT_POINT)
    
    tic;

    global SOLVER VARIABLE GEOMETRY
     
    %----------------------------------------------------------------------
    % Initial state, matrixes and load
    %----------------------------------------------------------------------    
    [BLCK,ste,ste_p,MAT_POINT,Disp_field,Int_var,Mat_state,GLOBAL,...
        stiff_mtx]=init(MAT_POINT);
    
    MATRIX=DYN_MATRIX; 
    
    [load_s,out_list1]=...
        calculate_forces(ste,MAT_POINT,Disp_field,Mat_state,MATRIX,BLCK);
    GLOBAL.OutputList(ste_p,:)=out_list1;
    
    %----------------------------------------------------------------------
    % Solver for each Block
    %---------------------------------------------------------------------- 
    while BLCK<SOLVER.BLOCKS+1
        
        if BLCK>1
            if SOLVER.Output(BLCK)~=SOLVER.Output(BLCK-1)  
                save(SOLVER.Output(BLCK), 'GEOMETRY', 'VARIABLE', 'SOLVER');
            end
        else
            save(SOLVER.Output(BLCK), 'GEOMETRY', 'VARIABLE', 'SOLVER');
        end

        if SOLVER.IMPLICIT(BLCK)==0
            [ste,ste_p,MAT_POINT,GLOBAL]=...
                Explicit_solver(BLCK,ste,ste_p,MAT_POINT,Disp_field,Int_var,...
                Mat_state,GLOBAL,MATRIX,load_s);
        else
            [ste,ste_p,MAT_POINT,GLOBAL]=...
                Implicit_solver(BLCK,ste,ste_p,MAT_POINT,Disp_field,Int_var,...
                Mat_state,GLOBAL,stiff_mtx,MATRIX,load_s);
        end
        BLCK=BLCK+1;
        
        if BLCK>SOLVER.BLOCKS
            break;
        else
            % Update
            [Disp_field,Mat_state,Int_var,stiff_mtx]=VECTORS.Update_ini(...
                BLCK,GLOBAL,ste,ste_p,Disp_field,Mat_state,Int_var,MAT_POINT);
            
            [load_s,out_list1]=...
                calculate_forces(ste,MAT_POINT,Disp_field,Mat_state,MATRIX,BLCK);
            GLOBAL.OutputList(ste_p,:)=out_list1;
        end
    end
    
    tfin=toc;

    disp(tfin);

end