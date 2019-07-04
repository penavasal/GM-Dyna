

function read_load

% File: read_load
%   Read loads from load.txt
%
% Date:
%   Version 1.0   18.04.2018

    global GEOMETRY SOLVER
    
    x_0=GEOMETRY.x_0;
    [~,sp]=size(x_0);
    
    %GEOM
    L=max(x_0(:,1));
    H=max(x_0(:,2));
    
    L0=min(x_0(:,1));
    H0=min(x_0(:,2));
    
    %FILE

    fid = fopen('load.txt', 'rt'); % opción rt para abrir en modo texto
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
        if strcmp(s1, '//')
                continue
        end
        if strcmp(s1,'LOAD')
            M=str2double(bb{t});
            if(M>=loads+1)
                break
            end
            s2=c{t};
            switch s2
                case 'POINT_LOAD'
                    TYPE(M)=1;
                    continue
                case 'LINE_LOAD'
                    TYPE(M)=2;
                    continue
                case 'VOLUME_ACCELERATION'
                    TYPE(M)=3;
                    continue
                case 'WATER_POINT_LOAD'
                    TYPE(M)=4;
                    continue
                case 'WATER_LINE_LOAD'
                    TYPE(M)=5;
                    continue
                otherwise
                    disp('Error, type of load not implented yet!')
                    stop
            end
        end
        if strcmp(s1,'X_RANGE')
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
        end
        if strcmp(s1,'Y_RANGE')
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
                elseif c{t}=='FULL'
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
        end
        if strcmp(s1,'VECTOR')
            VECTOR(1,M)=str2double(bb{t});
            VECTOR(2,M)=str2double(c{t});
            if sp==3
                VECTOR(3,M)=str2double(d{t});
            end
            continue
        end
        if strcmp(s1,'VALUE')
            VALUE(M)=bb{t};
            continue
        end
        if strcmp(s1,'INTERVAL')
            val=str2double(bb{t});
            if isnan(val)
                if bb{t}=='FULL' 
                    INTERVAL(1,M)=0;
                    INTERVAL(2,M)=SOLVER.Time_final;
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
        end
        
        fprintf('Error, unrecognized sequence: %s !! \n',s1)
        stop
        
    end
    
    fclose(fid);
    [load_mult_ini]=interval(INTERVAL,loads);
    [load_nds]=localization(RANGE,TYPE,loads);
    
    calculate_forces(load_mult_ini,load_nds,VALUE,VECTOR,TYPE,RANGE,loads);
end


function [Load_nds]=localization(R,T,mats)

    global GEOMETRY
    
    x_0=GEOMETRY.x_0;
    [nodes,~]=size(x_0);
    
    Load_nds=zeros(nodes,mats);
    for m=1:mats
        if T(m)==1 || T(m)==4
            x1=R(1,m*2-1);
            y1=R(2,m*2-1);
            d=zeros(nodes,1);
            for i=1:nodes
            	d(i)=sqrt((x_0(i,1)-x1)^2+(x_0(i,2)-y1)^2);
            end
            [~,j]=min(d);
            Load_nds(j,m)=1;   
        elseif T(m)==2 || T(m)==3 || T(m)==5
            for nodo=1:nodes
                tol=GEOMETRY.h_nds(nodo,1)/4;
                if (x_0(nodo,2)>=(R(2,m*2-1)-tol)) && (x_0(nodo,2)<=(R(2,m*2)+tol)) 
                    if (x_0(nodo,1)>=(R(1,m*2-1)-tol)) && (x_0(nodo,1)<=(R(1,m*2)+tol))
                        Load_nds(nodo,m)=1;
                    end
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

function calculate_forces(load_mult_ini,load_nds,VALUE,VECTOR,TYPE,RANGE,loads)
    
    global GEOMETRY VARIABLE SOLVER TIME LOAD
    
    g=VARIABLE.g;
    
    x_0=GEOMETRY.x_0;
    [nodes,sp]=size(x_0);
    
    LOAD.ext_forces_s = zeros(nodes*sp,loads);
    if SOLVER.UW==1
        LOAD.ext_forces_w = zeros(nodes*sp,loads);
    end
    LOAD.ext_acce = zeros(nodes*sp,loads);
    LOAD.load_mult = zeros(size(load_mult_ini));
    
    LOAD.size=loads;

    for m=1:loads
        
        dt=SOLVER.time_step;
        % Direction and nodes
        V=VECTOR(:,m)';
        nv=norm(V);
        if nv==0
            continue
        else
            V=V/norm(V);
        end
        
        % Values
        val=str2double(VALUE(m));
        if isnan(val)
            for i=1:SOLVER.step_final
                t=TIME.t(i);
                val=eval(VALUE(m));
                load_mult_ini(i,m)=load_mult_ini(i,m)*val;
            end
        else
            load_mult_ini(:,m)=load_mult_ini(:,m)*val;
        end
        if TYPE(m)==2 || TYPE(m)==5
            long=zeros(sp,1);
            for k=1:sp
                long(k)=(RANGE(k,m*2)-RANGE(k,m*2-1))^2;
            end
            area=sqrt(sum(long));
            if SOLVER.AXI
                if long(1) && long(2)==0
                    a1=RANGE(1,m*2);
                    a2=(RANGE(1,m*2)-area);
                    area=(a1^2-a2^2)*pi;
                elseif long(2) && long(1)==0
                    area=2*area*RANGE(1,m*2)*pi;
                end
            end
            LOAD.load_mult(:,m)=load_mult_ini(:,m)/area;
        else
            LOAD.load_mult(:,m)=load_mult_ini(:,m);
        end
        
        
        i=0;
        s=nnz(load_nds(:,m));
        nod_f=zeros(s,1);
        for nodo=1:nodes
            if load_nds(nodo,m)
                i=i+1;
                nod_f(i)=nodo;
            end
        end
        
        if TYPE(m)==1
            for j=1:i
                nn=nod_f(j);
                for k=1:sp
                    LOAD.ext_forces_s(nn*sp+1-k,m)=V(sp+1-k);
                end
            end
        elseif TYPE(m)==2
            if RANGE(1,m*2)-RANGE(1,m*2-1)==0
                [LOAD.ext_forces_s(:,m)]=dist_f(nod_f,x_0,V,SOLVER.AXI);
            elseif RANGE(2,m*2)-RANGE(2,m*2-1)==0
                [LOAD.ext_forces_s(:,m)]=dist_f_x(nod_f,x_0,V,SOLVER.AXI);
            else
                disp('not yet implemented, sorry');
            end
        elseif TYPE(m)==3
            for j=1:i
                nn=nod_f(j);
                for k=1:sp
                    LOAD.ext_acce(nn*sp+1-k,m)=V(sp+1-k);
                end
            end
        elseif TYPE(m)==4
            for j=1:i
                nn=nod_f(j);
                for k=1:sp
                    LOAD.ext_forces_w(nn*sp+1-k,m)=V(sp+1-k);
                end
            end
        elseif TYPE(m)==5
            if RANGE(1,m*2)-RANGE(1,m*2-1)==0
                [LOAD.ext_forces_w(:,m)]=dist_f(nod_f,x_0,V);
            elseif RANGE(2,m*2)-RANGE(2,m*2-1)==0
                [LOAD.ext_forces_w(:,m)]=dist_f_x(nod_f,x_0,V);
            else
                disp('not yet implemented, sorry');
            end
        end     
    end
end


function [ext_forces]=dist_f_x(nod,x_a,V,AXI)

    [nodes,sp]=size(x_a);
    [i,~]=size(nod);
    ext_forces=zeros(sp*nodes,1);

    maxy=1.0e-32;
    miny=1.0e32;
    for j=1:i
         if (x_a(nod(j),1)>maxy)
             maxy=x_a(nod(j),1);
         end
         if (x_a(nod(j),1)<miny)
             miny=x_a(nod(j),1);
         end
    end
    for j=1:i
        [d1,d2]=dist(x_a,nod,j,'X');
        if AXI
            rr=x_a(nod(j),1);
        end
        if (maxy==1.0e-32)
            d=0;
        elseif (x_a(nod(j),1)==maxy) || (x_a(nod(j),1)==miny)
            if AXI
                r=rr+d1;
                d=pi*r*abs(d1);
            else
                d=abs(d1)/2;
            end
        else
            if AXI
                r1=rr+d1;
                r2=rr+d2;
                d=pi*(r1*abs(d1)+r2*abs(d2));
            else
                d=abs(d1)/2+abs(d2)/2;
            end
        end
        f=V*d;
        ext_forces(nod(j)*sp-1)=f(1);
        ext_forces(nod(j)*sp)=f(2);
    end
end

function [ext_forces]=dist_f(nod,x_a,V,AXI)

    [nodes,sp]=size(x_a);
    [i,~]=size(nod);
    ext_forces=zeros(sp*nodes,1);

    maxy=1.0e-32;
    miny=1.0e32;
    for j=1:i
         if (x_a(nod(j),2)>maxy)
             maxy=x_a(nod(j),2);
         end
         if (x_a(nod(j),2)<miny)
             miny=x_a(nod(j),2);
         end
    end
    for j=1:i
        if AXI
            rr=x_a(nod(j),1);
        end
        [d1,d2]=dist(x_a,nod,j,'Y');
        if (maxy==1.0e-32)
            d=0;
        elseif (x_a(nod(j),2)==maxy) || (x_a(nod(j),2)==miny)
            d=d1/2;
        else
            d=d1/2+d2/2;
        end
        f=2*pi*rr*V*d;
        ext_forces(nod(j)*sp-1)=f(1);
        ext_forces(nod(j)*sp)=f(2);
    end
end

function [d1,d2]=dist(x_a,nd,j,r)
    x1=x_a(nd(j),1);
    y1=x_a(nd(j),2);
    d=zeros(length(nd),1);
    for i=1:length(nd)
        if(i~=j)
            d(i)=sqrt((x_a(nd(i),1)-x1)^2+(x_a(nd(i),2)-y1)^2);
            if strcmp(r,'X')
                if (x_a(nd(i),1)-x1)<0
                    d(i)=-d(i);
                end
            end
        else
            d(i)=1.0e32;
        end
    end
    d1=min(abs(d));
    t=0;
    i=0;
    while t==0&&i<length(nd)
        i=i+1;
        if abs(d(i))==abs(d1)
            t=1;
            d(i)=1.0e32;
        end  
    end
    d2=min(abs(d));   
end





