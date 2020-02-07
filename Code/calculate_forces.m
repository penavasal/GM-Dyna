function [load_s,OUTPUT]=...
            calculate_forces(STEP,MAT_POINT,Mat_state,MATRIX)

    global SOLVER LOAD GEOMETRY VARIABLE
    
    g=VARIABLE.g;
    t=STEP.t;
    BLCK=STEP.BLCK;

    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
    
    ext_forces_s = LOAD{BLCK}.ext_forces_s;
    ext_acce     = LOAD{BLCK}.ext_acce;
    value        = LOAD{BLCK}.value;
    if SOLVER.UW==1
        ext_forces_w = LOAD{BLCK}.ext_forces_w;
    end
    
    type=SOLVER.OutputType;
    [number,~]=size(type);
    OUTPUT(1,number)=0;
    
    load=zeros(GEOMETRY.nodes*df,LOAD{BLCK}.size);
    acce=zeros(GEOMETRY.nodes*df,LOAD{BLCK}.size);
    
    if LOAD{BLCK}.size>0
        if any(str2double(value(5,:)) == 3)
            [MATRIX]=MATRIX.lumped_mass_bf(MAT_POINT,Mat_state,MATRIX,BLCK);
        else
            MATRIX.l_mass=zeros(GEOMETRY.nodes*df);
        end
    end
    
    for m=1:LOAD{BLCK}.size
        
        if t<eval(value(2,m)) || t>eval(value(3,m))
            continue;
        end
        
        if strcmp(value(4,m),'VALUE')
            val=eval(value(1,m));
        elseif strcmp(value(4,m),'FUNCTION')
            val=eval(value(1,m));
        elseif strcmp(value(4,m),'FILE')
            [file_l]=load_file(value{1,m});
            val=interp1(file_l(:,1),file_l(:,2),[t],'pchip','extrap'); 
            if eval(value(5,m)) == 3
                val=val*g;
            end
        else
            disp('Error, wrong type of load!')
            stop
        end 
        vec=MATRIX.l_mass*ext_acce(:,m);
        
        for i=1:GEOMETRY.nodes 
            tt=1;
            if (SOLVER.UW==1 || SOLVER.UW==2) && SOLVER.IMPLICIT(BLCK)==1
                tt=-1;
            end
                
            if SOLVER.UW==1
            load((i-1)*df+1:(i-1)*df+sp,m)=load((i-1)*df+1:(i-1)*df+sp,m)...
                +tt*ext_forces_s((i-1)*sp+1:i*sp,m)*val;
            acce((i-1)*df+1:i*df,m)=acce((i-1)*df+1:i*df,m)...
                +tt*vec((i-1)*df+1:i*df,1)*val;
            load((i-1)*df+sp+1:i*df,m)=load((i-1)*df+sp+1:i*df,m)...
                +tt*ext_forces_w((i-1)*sp+1:i*sp,m)*val;
            elseif SOLVER.UW==0 || SOLVER.UW==2
            load((i-1)*df+1:(i-1)*df+sp,m)=load((i-1)*df+1:(i-1)*df+sp,m)...
                +tt*ext_forces_s((i-1)*sp+1:i*sp,m)*val;
            acce((i-1)*df+1:(i-1)*df+sp,m)=acce((i-1)*df+1:(i-1)*df+sp,m)...
                +tt*vec((i-1)*sp+1:i*sp,1)*val;
            end
        end
        
        for i=1:number
            if type(i,1)==2 && type(i,2)==m
                OUTPUT(1,i)=sum(load(:,m))+sum(acce(:,m));
            end
        end
    end
    
    load_s=sum(load,2)+sum(acce,2);
    
end

function[list]=load_file(ff)

    fid = fopen(ff, 'rt'); % opción rt para abrir en modo texto
    formato = '%s %s'; % formato de cada línea 
    data = textscan(fid, formato);
    
    list(:,1) = str2double(data{1});
    list(:,2) = str2double(data{2});

end
