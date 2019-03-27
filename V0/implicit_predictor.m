
function [Disp_field,Mat_state,Shape_function,FAIL]=implicit_predictor...
            (ste,GT,InvK,mass_mtx,damp_mtx,load_s,Disp_field,...
            Shape_function,Mat_state,Int_var,FAIL)
    
         % 1. Newton-Raphson
        [Disp_field,Mat_state,Shape_function,FAIL]=Newton_Raphson_solver...
            (ste,GT,InvK,mass_mtx,damp_mtx,load_s,Shape_function,...
            Disp_field,Int_var,Mat_state,FAIL);   
    
            % 2. Update variables
        [Disp_field.a,Disp_field.v]=G_solver_2...
            (Disp_field.d,Disp_field.a,Disp_field.v,ste);
        
        [Mat_state.xg]=update_mp(Disp_field.d,Shape_function,Mat_state.xg);
    
end