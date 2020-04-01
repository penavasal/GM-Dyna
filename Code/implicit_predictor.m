
function [Disp_field,Mat_state,MAT_POINT,STEP]=implicit_predictor...
            (STEP,stiff_mtx,mass_mtx,damp_mtx,load_s,Disp_field,...
            MAT_POINT,Mat_state,Int_var)
    
        global SOLVER
        
         % 1. Newton-Raphson
        [Disp_field,Mat_state,MAT_POINT,STEP]=Newton_Raphson_solver...
            (STEP,stiff_mtx,mass_mtx,damp_mtx,load_s,MAT_POINT,...
            Disp_field,Int_var,Mat_state);   
        
         % 2. Update variables
        if SOLVER.DYN(STEP.BLCK)==1  
            [Disp_field.a,Disp_field.v]=Time_Scheme.solver_2...
                (Disp_field.d,Disp_field.a,Disp_field.v,STEP);
        end
        
        [MAT_POINT]=update_mp(Disp_field.d,MAT_POINT);
    
end