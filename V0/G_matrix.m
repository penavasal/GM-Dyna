
function [matrix]=G_matrix(mass_mtx,stiff_mtx,damp_mtx,ste)

    global TI_param TIME

    af=TI_param(1);
    am=TI_param(2);
	delta=TI_param(3);
    alpha=TI_param(4);  
    theta=TI_param(5);
    
    if ste==1
        delta_t=TIME.t(ste+1)-TIME.t(ste);
    else
        delta_t=TIME.t(ste)-TIME.t(ste-1);
    end
    
    if alpha==0
        A=0;
        D=1/delta/delta_t;
    else
        A=(1-am)/alpha/delta_t/delta_t/theta/theta;
        D=(1-af)*delta/alpha/delta_t/theta;
    end

    
    matrix = (1-af)*stiff_mtx + A*mass_mtx + D*damp_mtx;
    
end


% delta=beta1=gamma
% alpha=beta2=beta
