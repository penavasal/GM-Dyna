
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
        
        save(SOLVER.Output(BLCK), 'GEOMETRY', 'VARIABLE', 'SOLVER');
        if SOLVER.DYN(BLCK)==0
            Static_solver(BLCK,ste,ste_p,MAT_POINT,Disp_field,Int_var,...
                    Mat_state,GLOBAL,stiff_mtx,MATRIX,load_s);
        else
            if SOLVER.IMPLICIT(BLCK)==0
                Explicit_solver(BLCK,ste,ste_p,MAT_POINT,Disp_field,Int_var,...
                    Mat_state,GLOBAL,MATRIX,load_s)
            else
                Implicit_solver(BLCK,ste,ste_p,MAT_POINT,Disp_field,Int_var,...
                    Mat_state,GLOBAL,stiff_mtx,MATRIX,load_s)
            end
        end
        BLCK=BLCK+1;
    end


    tfin=toc;

    disp(tfin);

end