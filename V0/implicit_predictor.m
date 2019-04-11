
function [Disp_field,Mat_state,MAT_POINT,FAIL]=implicit_predictor...
            (ste,GT,InvK,mass_mtx,damp_mtx,load_s,Disp_field,...
            MAT_POINT,Mat_state,Int_var,FAIL)
    
         % 1. Newton-Raphson
        [Disp_field,Mat_state,MAT_POINT,FAIL]=Newton_Raphson_solver...
            (ste,GT,InvK,mass_mtx,damp_mtx,load_s,MAT_POINT,...
            Disp_field,Int_var,Mat_state,FAIL);   
    
            % 2. Update variables
        [Disp_field.a,Disp_field.v]=G_solver_2...
            (Disp_field.d,Disp_field.a,Disp_field.v,ste);
        
        [MAT_POINT]=update_mp(Disp_field.d,MAT_POINT);
    
end