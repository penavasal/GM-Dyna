function [load_s,OUTPUT]=...
            calculate_forces(ste,MAT_POINT,Disp_field,Mat_state,OUTPUT,MATRIX)

    global SOLVER LOAD GEOMETRY

    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
    
    ext_forces_s = LOAD.ext_forces_s;
    ext_acce     = LOAD.ext_acce;
    load_mult    = LOAD.load_mult;
    if SOLVER.UW==1
        ext_forces_w = LOAD.ext_forces_w;
    end

    
    [MATRIX]=MATRIX.lumped_mass_bf(MAT_POINT,Mat_state,MATRIX);
    
    load=zeros(GEOMETRY.nodes*df,LOAD.size);
    acce=zeros(GEOMETRY.nodes*df,LOAD.size);
    
    for m=1:LOAD.size
        vec=MATRIX.l_mass*ext_acce(:,m);
        for i=1:GEOMETRY.nodes 
%             if SOLVER.AXI
%                 %if Disp_field.x_a(i,1)==0
%                 %    t=2*pi*GEOMETRY.h_nds(i)/20;  
%                 %else
%                     t=2*pi*Disp_field.x_a(i,1);
%                 %end
%             else
%                 t=1;
%             end
            t=1;
            if SOLVER.UW==1 && (SOLVER.IMPLICIT==1 || SOLVER.step0==1) 
                t=-1;
            end
            
            if SOLVER.UW==2 && (SOLVER.IMPLICIT==1 || SOLVER.step0==1) 
                t=-1;
            end
                
            if SOLVER.UW==1
            load((i-1)*df+1:(i-1)*df+sp,m)=load((i-1)*df+1:(i-1)*df+sp,m)...
                +t*ext_forces_s((i-1)*sp+1:i*sp,m)*load_mult(ste,m);
            acce((i-1)*df+1:i*df,m)=acce((i-1)*df+1:i*df,m)...
                +t*vec((i-1)*df+1:i*df,1)*load_mult(ste,m);
            load((i-1)*df+sp+1:i*df,m)=load((i-1)*df+sp+1:i*df,m)...
                +t*ext_forces_w((i-1)*sp+1:i*sp,m)*load_mult(ste,m);
            elseif SOLVER.UW==0 || SOLVER.UW==2
            load((i-1)*df+1:(i-1)*df+sp,m)=load((i-1)*df+1:(i-1)*df+sp,m)...
                +t*ext_forces_s((i-1)*sp+1:i*sp,m)*load_mult(ste,m);
            acce((i-1)*df+1:(i-1)*df+sp,m)=acce((i-1)*df+1:(i-1)*df+sp,m)...
                +t*vec((i-1)*sp+1:i*sp,1)*load_mult(ste,m);
            end
        end
        
        for i=1:OUTPUT.number
            if OUTPUT.type(i,1)==1 && OUTPUT.type(i,2)==m
                OUTPUT.inst(i)=sum(load(:,m));
            end
        end
    end
    
    load_s=sum(load,2)+sum(acce,2);
    
end
