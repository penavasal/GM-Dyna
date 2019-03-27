function TIME_INT_var

    global MATERIAL VARIABLE SOLVER GEOMETRY TIME

    
    h=GEOMETRY.h_ini;
    tt(GEOMETRY.elements,1)=0;
    for e=1:GEOMETRY.elements
        if SOLVER.UW
            tt(e,1)=min(h(e)/MATERIAL.MAT(6,MATERIAL.e(e)),...
            h(e)/sqrt(MATERIAL.MAT(28,MATERIAL.e(e))/VARIABLE.rho_w));
        else
            tt(e,1)=h(e)/MATERIAL.MAT(6,MATERIAL.e(e));
        end
    end
    TT=min(tt);
    CFL=SOLVER.time_step/TT;
    fprintf('%f of CFL\n',CFL);

    delta_t=SOLVER.time_step;
    
    ste=1;
    t(ste,1)=0;
    while t(ste,1)<SOLVER.Time_final
        ste=ste+1;
        delta_t=delta_t*SOLVER.time_factor;
        t(ste,1)=t(ste-1,1)+delta_t;
    end
    SOLVER.Time_final=t(ste,1);
    SOLVER.step_final=ste;
    TIME.t=t;
    

    if SOLVER.SAVE_I==1
        SOLVER.dim=floor(SOLVER.step_final/SOLVER.SAVE_I);
    else
        SOLVER.dim=floor(SOLVER.step_final/SOLVER.SAVE_I)+1;
    end
    TIME.tp=zeros(SOLVER.dim,1);
    fprintf('%i plot steps\n',SOLVER.dim);
    fprintf('Save %i times before the final\n',round(SOLVER.dim/SOLVER.SAVE_F));

end