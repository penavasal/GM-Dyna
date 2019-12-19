
function [Disp_field,Mat_state,MAT_POINT]=implicit_predictor...
            (ste,stiff_mtx,mass_mtx,damp_mtx,load_s,Disp_field,...
            MAT_POINT,Mat_state,Int_var,BLCK)
    
         % 1. Newton-Raphson
        [Disp_field,Mat_state,MAT_POINT]=Newton_Raphson_solver...
            (ste,stiff_mtx,mass_mtx,damp_mtx,load_s,MAT_POINT,...
            Disp_field,Int_var,Mat_state,BLCK);   
    
            % 2. Update variables
        [Disp_field.a,Disp_field.v]=Time_Scheme.solver_2...
            (Disp_field.d,Disp_field.a,Disp_field.v,ste,BLCK);
        
        [MAT_POINT]=update_mp(Disp_field.d,MAT_POINT);
    
end