
function read_boundary

% File: read_boundary
%   Read loads from boundary.txt
%
% Date:
%   Version 1.0   19.04.2018

    global GEOMETRY SOLVER
    
    x_0=GEOMETRY.x_0;
    [~,sp]=size(x_0);
    
    %GEOM
    L=max(x_0(:,1));
    H=max(x_0(:,2));
    
    L0=min(x_0(:,1));
    H0=min(x_0(:,2));
    
    %FILE

    fid = fopen('boundary.txt', 'rt'); % opción rt para abrir en modo texto
    formato = '%s %s %s %s'; % formato de cada línea 
    data = textscan(fid, formato, 'HeaderLines', 1);
    
    a = data{1};
    % Convertir a vector numérico
    [long,~]=size(a);
    bb = data{2};
    c = data{3};
    d = data{4};
    
    l=1;
    loads=str2double(bb{l});
    while isnan(loads)
        l=l+1;
        loads=str2double(bb{l});
    end
    RANGE=zeros(sp,2*loads); % Range of loads
    VECTOR=zeros(sp,loads); % Directions of loads
    VALUE = strings(loads);
    INTERVAL=zeros(2,loads); % Intervals of loads
    TYPE=zeros(loads,1);

    M=0;
    t=l;
    
    while (M<loads+1) && (t<long) 
        
        t=t+1;
        s1=a{t};
        switch s1
            case '//'
                continue
            case 'BOUNDARY'
                M=str2double(bb{t});
                if(M>=loads+1)
                    break
                end
                s2=c{t};
                switch s2
                    case 'DISPLACEMENT'
                        TYPE(M)=1;
                        continue
                    case 'WATER_DISPLACEMENT'
                        TYPE(M)=2;
                        continue
                    case 'VELOCITY'
                        TYPE(M)=3;
                        continue
                    case 'WATER_VELOCITY'
                        TYPE(M)=4;
                        continue
                    otherwise
                        disp('Error, type of boundary not implemented yet!')
                        stop
                end
            case 'X_RANGE'
                cc = str2double(bb{t});
                if isnan(cc)
                    if strcmp(bb{t},'FULL')
                        RANGE(1,M*2-1)=L;
                    elseif strcmp(bb{t},'INI')
                        RANGE(1,M*2-1)=L0;
                    else
                        disp('Error, wrong load X range!')
                    end
                else
                    RANGE(1,M*2-1)=cc;
                end
                cc = str2double(c{t});
                if isnan(cc)
                    len=strlength(c{t});
                    if len==0 
                        RANGE(1,M*2)=RANGE(1,M*2-1);
                    elseif strcmp(c{t},'FULL')
                        RANGE(1,M*2)=L;
                    elseif strcmp(c{t},'INI')
                        RANGE(1,M*2)=L0;
                    else
                        disp('Error, wrong load X range!')
                        stop
                    end
                else
                    RANGE(1,M*2)=cc;
                end
                continue
            case 'Y_RANGE'
                cc = str2double(bb{t});
                if isnan(cc)
                    if strcmp(bb{t},'FULL')
                        RANGE(2,M*2-1)=H;
                    elseif strcmp(bb{t},'INI')
                        RANGE(2,M*2-1)=H0;
                    else
                        disp('Error, wrong load Y range!')
                    end
                else
                    RANGE(2,M*2-1)=cc;
                end
                cc = str2double(c{t});
                if isnan(cc)
                    len=strlength(c{t});
                    if len==0 
                        RANGE(2,M*2)=RANGE(2,M*2-1);
                    elseif strcmp(c{t},'FULL')
                        RANGE(2,M*2)=H;
                    elseif strcmp(c{t},'INI')
                        RANGE(2,M*2)=H0;
                    else
                        disp('Error, wrong load Y range!')
                        stop
                    end
                else
                    RANGE(2,M*2)=cc;
                end
                continue
            case 'VECTOR'
                VECTOR(1,M)=str2double(bb{t});
                VECTOR(2,M)=str2double(c{t});
                if sp==3
                    VECTOR(3,M)=str2double(d{t});
                end
                continue
            case 'VALUE'
                VALUE(M)=bb{t};
                continue
            case 'INTERVAL'
                val=str2double(bb{t});
                val2=str2double(c{t});
                if isnan(val)
                    if strcmp(bb{t},'FULL')
                        INTERVAL(1,M)=0;
                        INTERVAL(2,M)=SOLVER.Time_final;
                    elseif strcmp(bb{t},'INI')
                        INTERVAL(1,M)=0;
                        if isnan(val2)
                            if strcmp(c{t},'FULL')
                                INTERVAL(2,M)=SOLVER.Time_final;
                            elseif strcmp(c{t},'INI') || strcmp(c{t},'')
                                INTERVAL(2,M)=0;
                            end
                        else
                            if val2==0
                                INTERVAL(2,M)=0;
                            else
                                INTERVAL(2,M)=val2;
                            end
                        end
                    else
                        disp('Error, wrong load interval!')
                        stop
                    end
                else
                    INTERVAL(1,M)=val;
                    val=str2double(c{t});
                    if isnan(val)
                        if c{t}=='FULL'
                            INTERVAL(2,M)=SOLVER.Time_final;
                        else
                            disp('Error, wrong load interval!')
                            stop
                        end
                    else
                        INTERVAL(2,M)=val;
                    end
                end
                continue
            otherwise
            	fprintf('Error, unrecognized sequence: %s !! \n',s1)
                stop
        end
                
    end
    
    fclose(fid);
    
    [b_mult_ini]=interval(INTERVAL,loads);
    [b_nds]=localization(RANGE,loads);
    
    calculate_boundaries(b_mult_ini,b_nds,VALUE,VECTOR,TYPE,loads);
end

function [Load_nds]=localization(R,mats)

    global GEOMETRY
    x_0=GEOMETRY.x_0;
    [nodes,~]=size(x_0);
    Load_nds=zeros(nodes,mats);
    for m=1:mats
        for nodo=1:nodes
            tol=GEOMETRY.h_nds(nodo,1)/5;
            if (x_0(nodo,2)>=R(2,m*2-1)-tol) && (x_0(nodo,2)<=R(2,m*2)+tol) 
                if (x_0(nodo,1)>=R(1,m*2-1)-tol) && (x_0(nodo,1)<=R(1,m*2)+tol)
                    Load_nds(nodo,m)=1;
                end
            end
        end
    end

end

function [load_mult]=interval(INTERVAL,loads)

    global SOLVER TIME

    load_mult=zeros(SOLVER.step_final,loads);
    
    for m=1:loads
        for i=1:SOLVER.step_final
            t=TIME.t(i);
            if t>=INTERVAL(1,m) && t<=INTERVAL(2,m)
                load_mult(i,m)=1;
            end
        end
        
    end

end

function calculate_boundaries(load_mult_ini,load_nds,VALUE,VECTOR,TYPE,loads)
    
    global GEOMETRY SOLVER TIME BOUNDARY
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    
    BOUNDARY.constrains = zeros(GEOMETRY.nodes*df,loads);
    BOUNDARY.dad = zeros(GEOMETRY.nodes*df,loads);
    BOUNDARY.vad = zeros(GEOMETRY.nodes*df,loads);
    b_mult = strings(size(load_mult_ini)); 
    
    BOUNDARY.size=loads;

    for m=1:loads
        
        % Direction
        V=VECTOR(:,m)';
        nv=norm(V);
        if nv==0
            continue
        else
            V=V/norm(V);
        end
        
        % Nodes
        for i=1:GEOMETRY.nodes
            if load_nds(i,m)
                if SOLVER.UW==0 && (TYPE(m)==2 || TYPE(m)==4)
                    disp('error, take care of the water boundary conditions!!');
                elseif (TYPE(m)==2) || (TYPE(m)==1 && SOLVER.UW==0)
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY.dad(i*df+1-k,m)=V(sp+1-k);
                            BOUNDARY.constrains(i*df+1-k,m)=1;
                        end
                    end
                elseif TYPE(m)==1 && SOLVER.UW==1
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY.dad(i*df-sp+1-k,m)=V(sp+1-k);
                            BOUNDARY.constrains(i*df-sp+1-k,m)=1;
                        end
                    end
                elseif (TYPE(m)==4) || (TYPE(m)==3 && SOLVER.UW==0)
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY.vad(i*df+1-k,m)=V(sp+1-k);
                            BOUNDARY.constrains(i*df+1-k,m)=2;
                        end
                    end
                elseif TYPE(m)==3 && SOLVER.UW==1
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY.vad(i*df-sp+1-k,m)=V(sp+1-k);
                            BOUNDARY.constrains(i*df-sp+1-k,m)=2;
                        end
                    end
                end
                if SOLVER.UW==0 && sp~=df
                    % No water
                    BOUNDARY.constrains(i*df,m)=1;
                    BOUNDARY.constrains(i*df-1,m)=1;
                end
            end
        end
        
        %Values
        for i=1:SOLVER.step_final
            val=str2double(VALUE(m));
            if isnan(val)
                t=TIME.t(i);
                val=eval(VALUE(m));
            end
            if load_mult_ini(i,m)==1
                b_mult(i,m)=num2str(val);
            else
                b_mult(i,m)='NULL';
            end
        end
    end
    BOUNDARY.b_mult=b_mult;
end