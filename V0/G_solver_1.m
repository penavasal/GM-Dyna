function [incr_d]=G_solver_1(InvK,FT,a1,v1,ste)

    global TI_param TIME

    %af=TI_param(1);
    %am=TI_param(2);
	%delta=TI_param(3);
    alpha=TI_param(4);
    theta=TI_param(5);
    
    if ste==1
        delta_t=TIME.t(ste+1)-TIME.t(ste);
    else
        delta_t=TIME.t(ste)-TIME.t(ste-1);
    end
    
    if theta==1
        incr_d=InvK*FT;
    else
        A=1/alpha/delta_t/delta_t/theta/theta;
        B=1/alpha/delta_t/theta;
        C=1/2/alpha;
        
        incr_d_th=InvK*FT;
        incr_a_th=A*incr_d_th-B*v1(:,2)-C*a1(:,2);

        incr_d=incr_d_th-(delta_t*(theta-1)*v1(:,2)+delta_t^2/2*(theta^2-1)*a1(:,2)...
                          +alpha*delta_t^2*(theta^2-1)*incr_a_th);
    end
    
%     a1(:,1)=a1(:,2)+incr_a_th/theta;
%     v1(:,1)=v1(:,2)+delta_t*a1(:,2)+delta*delta_t/theta*incr_a_th;
%     d1(:,1)=d1(:,2)+delta_t*v1(:,2)+delta_t^2/2*a1(:,2)+...
%         alpha*delta_t^2/theta*incr_a_th;
    
end