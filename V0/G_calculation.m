function [GT]=G_calculation(d1,a1,v1,Fint,Mm,Cm,load,load1,ste)
   
    
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
        B=0;
        C=0;
        D=1/delta/delta_t;
        E=1-1/delta;
        F=0;
    else
        A=(1-am)/alpha/delta_t/delta_t/theta/theta;
        B=(1-am)/alpha/delta_t/theta;
        C=(1-am)/2/alpha-1;
        D=(1-af)*delta/alpha/delta_t/theta;
        E=1-(1-af)*delta/alpha;
        F=(1-af)*(1-delta/2/alpha)*delta_t*theta;
    end
    
    G=theta*(1-af);
    H=(1-af);
   
    du=d1(:,1)-d1(:,2);

    GT= G*(load-load1)+load1 - H*Fint(:,1)...
            -Mm*(A*du-B*v1(:,2)-C*a1(:,2))...
            -Cm*(D*du+E*v1(:,2)+F*a1(:,2));

end

