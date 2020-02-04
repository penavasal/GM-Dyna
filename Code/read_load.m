

function read_load(filetxt,BLCK,NODE_LIST)

% File: read_load
%   Read loads from load.txt & .dat
%
% Date:
%   Version 2.0   25.11.2019

    global GEOMETRY SOLVER
    
    x_0=GEOMETRY.x_0;
    [~,sp]=size(x_0);
    
    %GEOM
    L=max(x_0(:,1));
    H=max(x_0(:,2));
    
    L0=min(x_0(:,1));
    H0=min(x_0(:,2));
    
    %FILE

    fid = fopen(filetxt, 'rt'); % opción rt para abrir en modo texto
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
        
    %RANGE=zeros(sp,2*loads); % Range of loads
    NLIST=strings(loads,1); % Lista de nodos
    VECTOR=zeros(sp,loads); % Directions of loads
    VALUE = strings(loads,1);
    OUT = zeros(loads,1);
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
%         if strcmp(s1,'X_RANGE')
%             cc = str2double(bb{t});
%             if isnan(cc)
%                 if strcmp(bb{t},'FULL')
%                     RANGE(1,M*2-1)=L;
%                 elseif strcmp(bb{t},'INI')
%                     RANGE(1,M*2-1)=L0;
%                 else
%                     disp('Error, wrong load X range!')
%                 end
%             else
%                 RANGE(1,M*2-1)=cc;
%             end
%             cc = str2double(c{t});
%             if isnan(cc)
%                 len=strlength(c{t});
%                 if len==0 
%                     RANGE(1,M*2)=RANGE(1,M*2-1);
%                 elseif strcmp(c{t},'FULL')
%                     RANGE(1,M*2)=L;
%                 elseif strcmp(c{t},'INI')
%                     RANGE(1,M*2)=L0;
%                 else
%                     disp('Error, wrong load X range!')
%                     stop
%                 end
%             else
%                 RANGE(1,M*2)=cc;
%             end
%             continue
%         end
%         if strcmp(s1,'Y_RANGE')
%             cc = str2double(bb{t});
%             if isnan(cc)
%                 if strcmp(bb{t},'FULL')
%                     RANGE(2,M*2-1)=H;
%                 elseif strcmp(bb{t},'INI')
%                     RANGE(2,M*2-1)=H0;
%                 else
%                     disp('Error, wrong load Y range!')
%                 end
%             else
%                 RANGE(2,M*2-1)=cc;
%             end
%             cc = str2double(c{t});
%             if isnan(cc)
%                 len=strlength(c{t});
%                 if len==0 
%                     RANGE(2,M*2)=RANGE(2,M*2-1);
%                 elseif c{t}=='FULL'
%                     RANGE(2,M*2)=H;
%                 elseif strcmp(c{t},'INI')
%                 	RANGE(2,M*2)=H0;
%                 else
%                     disp('Error, wrong load Y range!')
%                     stop
%                 end
%             else
%                 RANGE(2,M*2)=cc;
%             end
%             continue
%         end
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
        if strcmp(s1,'OUTPUT')
            OUT(M)=str2double(bb{t});
            continue
        end
        if strcmp(s1,'NODE_LIST')
            NLIST(M)=bb{t};
            continue
        end
        if strcmp(s1,'INTERVAL')
            val=str2double(bb{t});
            if isnan(val)
                if bb{t}=='FULL' 
                    if BLCK==1
                        INTERVAL(1,M)=0;
                    else
                        INTERVAL(1,M)=SOLVER.Time_final(BLCK-1);
                    end
                    INTERVAL(2,M)=SOLVER.Time_final(BLCK);
                else
                    disp('Error, wrong load interval!')
                    stop
                end
            else
                INTERVAL(1,M)=val;
                val=str2double(c{t});
                if isnan(val)
                    if c{t}=='FULL'
                        INTERVAL(2,M)=SOLVER.Time_final(BLCK);
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
    
    interval(INTERVAL,loads,BLCK,VALUE,TYPE);
    %[load_nds]=localization(RANGE,TYPE,loads);
    
    calculate_forces(NODE_LIST,VECTOR,TYPE,NLIST,loads,BLCK);
    
    for i=1:loads
        if OUT(i)==1
            tt=0;
            for j=1:length(SOLVER.OutputType(:,1))
                if SOLVER.OutputType(j,1)==0
                    SOLVER.OutputType(j,1)=2;
                    SOLVER.OutputType(j,2)=i;
                    tt=0;
                    break;
                elseif SOLVER.OutputType(j,1)==2 && SOLVER.OutputType(j,2)==i
                    tt=0;
                    break;
                else
                    tt=1;
                end
            end
            if tt==1
                j=j+1;
                SOLVER.OutputType(j,1)=2;
                SOLVER.OutputType(j,2)=i;
            end
        end
    end
    
end

function interval(INTERVAL,loads,BLCK,VALUE,TYPE)

    global LOAD
        
    %LOAD{BLCK}.value(5,loads)=0;
    LOAD{BLCK}.value=strings(5,loads);
    
    for m=1:loads
        % Value
        val=str2double(VALUE(m));
        if isnan(val)
            if isfile(strcat(VALUE(m),'.txt'))
                ff=strcat(VALUE(m),'.txt');
                LOAD{BLCK}.value(1,m)={load_file(ff)};
                LOAD{BLCK}.value(4,m)='FILE';
            else
                LOAD{BLCK}.value(4,m)='FUNCTION';
                LOAD{BLCK}.value(1,m)=VALUE(m);
            end
        else
            LOAD{BLCK}.value(4,m)='VALUE';
            LOAD{BLCK}.value(1,m)=val;
        end
        LOAD{BLCK}.value(2,m)=INTERVAL(1,m);
        LOAD{BLCK}.value(3,m)=INTERVAL(2,m);
        LOAD{BLCK}.value(5,m)=TYPE(m);
    end
end

function calculate_forces...
    (NODE_LIST,VECTOR,TYPE,NLIST,loads,BLCK)
    
    global GEOMETRY SOLVER LOAD
    
       
    x_0=GEOMETRY.x_0;
    [nodes,sp]=size(x_0);
    df=GEOMETRY.df;
    
    if loads
        LOAD{BLCK}.ext_forces_s = zeros(nodes*sp,loads);
        if SOLVER.UW==1
            LOAD{BLCK}.ext_forces_w = zeros(nodes*sp,loads);
        end
        LOAD{BLCK}.ext_acce = zeros(nodes*df,loads);    
        LOAD{BLCK}.size=loads;

        pls=NODE_LIST.pls;
        lls=NODE_LIST.lls;
        vls=NODE_LIST.vls;

        PL=NODE_LIST.PL;
        LL=NODE_LIST.LL;
        VL=NODE_LIST.VL;  
    else
        LOAD{BLCK}.ext_forces_s = 0;
        if SOLVER.UW==1
            LOAD{BLCK}.ext_forces_w = 0;
        end
        LOAD{BLCK}.ext_acce = 0;   
        LOAD{BLCK}.size=loads;
    end

    for m=1:loads
        
        % Direction and nodes
        V=VECTOR(:,m)';
        nv=norm(V);
        if nv==0
            continue
        else
            V=V/norm(V);
        end
                
        %Area
%         if TYPE(m)==2 || TYPE(m)==5
%             long=zeros(sp,1);
%             for k=1:sp
%                 long(k)=(RANGE(k,m*2)-RANGE(k,m*2-1))^2;
%             end
%             area=sqrt(sum(long));
%             if SOLVER.AXI
%                 if long(1) && long(2)==0
%                     a1=RANGE(1,m*2);
%                     a2=(RANGE(1,m*2)-area);
%                     area=(a1^2-a2^2)*pi;
%                 elseif long(2) && long(1)==0
%                     area=2*area*RANGE(1,m*2)*pi;
%                 end
%             end
%             LOAD.load_mult(:,m)=load_mult_ini(:,m)/area;
%         else
%             LOAD.load_mult(:,m)=load_mult_ini(:,m);
%         end
        
        
%         i=0;
%         s=nnz(load_nds(:,m));
%         nod_f=zeros(s,1);
%         for nodo=1:nodes
%             if load_nds(nodo,m)
%                 i=i+1;
%                 nod_f(i)=nodo;
%             end
%         end


        %NODE LIST
        
        ll=str2double(NLIST(m));
        if isnan(ll)
            if strcmp(NLIST(m),'FULL') && (TYPE(m)~=2 || TYPE(m)~=5)
                nod_f=linspace(1,nodes);
            else
                fprintf('Error, unrecognized list of nodes: %s !! \n',ll)
                stop
            end
        else
            if TYPE(m)==1 || TYPE(m)==4
                if ll>pls
                    fprintf('Error, unrecognized list of nodes: %i !! \n',ll)
                    stop
                else
                    nod_f=PL{ll};
                end 
            elseif TYPE(m)==3
                if ll>vls
                    fprintf('Error, unrecognized list of nodes: %i !! \n',ll)
                    stop
                else
                    nod_f=VL{ll};
                end
            elseif TYPE(m)==2 || TYPE(m)==5
                if ll>lls
                    fprintf('Error, unrecognized list of nodes: %i !! \n',ll)
                    stop
                else
                    nod_f=LL{ll};
                end
            end
        end

        % LOADING
        [long,~]=size(nod_f);
        for j=1:long
            if TYPE(m)==1
                nn=nod_f(j);
                for k=1:sp
                    LOAD{BLCK}.ext_forces_s(nn*sp+1-k,m)=V(sp+1-k);
                end
            elseif TYPE(m)==2
                if SOLVER.AXI
                    rr=(x_0(nod_f(j,1),1)+x_0(nod_f(j,2),1))/2;
                end
                area=zeros(sp,1);
                for k=1:sp
                    area(k)=(x_0(nod_f(j,1),k)-x_0(nod_f(j,2),k))^2;
                end
                d=sqrt(sum(area));
                if SOLVER.AXI
                    f=2*pi*rr*V*d/2;
                else
                    f=V*d/2;
                end
                LOAD{BLCK}.ext_forces_s(nod_f(j,1)*sp-1,m)=...
                    LOAD{BLCK}.ext_forces_s(nod_f(j,1)*sp-1,m)+f(1);
                LOAD{BLCK}.ext_forces_s(nod_f(j,1)*sp,m)=...
                    LOAD{BLCK}.ext_forces_s(nod_f(j,1)*sp,m)+f(2);
                LOAD{BLCK}.ext_forces_s(nod_f(j,2)*sp-1,m)=...
                    LOAD{BLCK}.ext_forces_s(nod_f(j,2)*sp-1,m)+f(1);
                LOAD{BLCK}.ext_forces_s(nod_f(j,2)*sp,m)=...
                    LOAD{BLCK}.ext_forces_s(nod_f(j,2)*sp,m)+f(2);
%             if RANGE(1,m*2)-RANGE(1,m*2-1)==0
%                 [LOAD.ext_forces_s(:,m)]=dist_f(nod_f,x_0,V,SOLVER.AXI);
%             elseif RANGE(2,m*2)-RANGE(2,m*2-1)==0
%                 [LOAD.ext_forces_s(:,m)]=dist_f_x(nod_f,x_0,V,SOLVER.AXI);
%             else
%                 disp('not yet implemented, sorry');
%             end
            elseif TYPE(m)==3
                nn=nod_f(j);
                if SOLVER.UW==0
                    for k=1:sp
                        LOAD{BLCK}.ext_acce(nn*sp+1-k,m)=V(sp+1-k);
                    end
                elseif SOLVER.UW==1
                    for k=1:sp
                        LOAD{BLCK}.ext_acce(nn*df+1-k,m)=V(sp+1-k);
                        LOAD{BLCK}.ext_acce(nn*df-sp+1-k,m)=V(sp+1-k);
                    end
                elseif SOLVER.UW==2
                    for k=1:sp
                        LOAD{BLCK}.ext_acce(nn*df-k,m)=V(sp+1-k);
                    end
                end
            elseif TYPE(m)==4
                nn=nod_f(j);
                for k=1:sp
                    LOAD{BLCK}.ext_forces_w(nn*sp+1-k,m)=V(sp+1-k);
                end
            elseif TYPE(m)==5
                if SOLVER.AXI
                    rr=(x_0(nod_f(j,1),1)+x_0(nod_f(j,2),1))/2;
                end
                area=zeros(sp,1);
                for k=1:sp
                    area(k)=(x_0(nod_f(j,1),1)+x_0(nod_f(j,2),1))^2;
                end
                d=sqrt(sum(area));
                if SOLVER.AXI
                    f=2*pi*rr*V*d/2;
                else
                    f=V*d/2;
                end
                LOAD{BLCK}.ext_forces_s(nod_f(j,1)*sp-1,m)=...
                    LOAD{BLCK}.ext_forces_w(nod_f(j,1)*sp-1,m)+f(1);
                LOAD{BLCK}.ext_forces_s(nod_f(j,1)*sp,m)=...
                    LOAD{BLCK}.ext_forces_w(nod_f(j,1)*sp,m)+f(2);
                LOAD{BLCK}.ext_forces_s(nod_f(j,2)*sp-1,m)=...
                    LOAD{BLCK}.ext_forces_w(nod_f(j,2)*sp-1,m)+f(1);
                LOAD{BLCK}.ext_forces_s(nod_f(j,2)*sp,m)=...
                    LOAD{BLCK}.ext_forces_w(nod_f(j,2)*sp,m)+f(2);
            end
        end     
    end
end

function[list]=load_file(ff)

    fid = fopen(ff, 'rt'); % opción rt para abrir en modo texto
    formato = '%s %s'; % formato de cada línea 
    data = textscan(fid, formato, 'HeaderLines', 1);
    
    list(:,1) = str2double(data{1});
    list(:,2) = str2double(data{2});

end

% 
% function [d1,d2]=dist(x_a,nd,j,r)
%     x1=x_a(nd(j),1);
%     y1=x_a(nd(j),2);
%     d=zeros(length(nd),1);
%     for i=1:length(nd)
%         if(i~=j)
%             d(i)=sqrt((x_a(nd(i),1)-x1)^2+(x_a(nd(i),2)-y1)^2);
%             if strcmp(r,'X')
%                 if (x_a(nd(i),1)-x1)<0
%                     d(i)=-d(i);
%                 end
%             end
%         else
%             d(i)=1.0e32;
%         end
%     end
%     [~,i]=min(abs(d));
%     d1=d(i);
%     d(i)=1.0e32;
%     d2=min(abs(d));   
% end

% 
% function [ext_forces]=dist_f_x(nod,x_a,V,AXI)
% 
%     [nodes,sp]=size(x_a);
%     [i,~]=size(nod);
%     ext_forces=zeros(sp*nodes,1);
% 
%     maxy=1.0e-32;
%     miny=1.0e32;
%     for j=1:i
%          if (x_a(nod(j),1)>maxy)
%              maxy=x_a(nod(j),1);
%          end
%          if (x_a(nod(j),1)<miny)
%              miny=x_a(nod(j),1);
%          end
%     end
%     for j=1:i
%         [d1,d2]=dist(x_a,nod,j,'X');
%         if AXI
%             rr=x_a(nod(j),1);
%         end
%         if (maxy==1.0e-32)
%             d=0;
%         elseif (x_a(nod(j),1)==maxy) || (x_a(nod(j),1)==miny)
%             if AXI
%                 r=rr+d1/2;
%                 d=pi*r*abs(d1);
%             else
%                 d=abs(d1)/2;
%             end
%         else
%             if AXI
%                 r1=rr+d1/2;
%                 r2=rr+d2/2;
%                 d=pi*(r1*abs(d1)+r2*abs(d2));
%             else
%                 d=abs(d1)/2+abs(d2)/2;
%             end
%         end
%         f=V*d;
%         ext_forces(nod(j)*sp-1)=f(1);
%         ext_forces(nod(j)*sp)=f(2);
%     end
% end
% 
% function [ext_forces]=dist_f(nod,x_a,V,AXI)
% 
%     [nodes,sp]=size(x_a);
%     [i,~]=size(nod);
%     ext_forces=zeros(sp*nodes,1);
% 
%     maxy=1.0e-32;
%     miny=1.0e32;
%     for j=1:i
%          if (x_a(nod(j),2)>maxy)
%              maxy=x_a(nod(j),2);
%          end
%          if (x_a(nod(j),2)<miny)
%              miny=x_a(nod(j),2);
%          end
%     end
%     for j=1:i
%         if AXI
%             rr=x_a(nod(j),1);
%         end
%         [d1,d2]=dist(x_a,nod,j,'Y');
%         if (maxy==1.0e-32)
%             d=0;
%         elseif (x_a(nod(j),2)==maxy) || (x_a(nod(j),2)==miny)
%             d=d1/2;
%         else
%             d=d1/2+d2/2;
%         end
%         if AXI
%             f=2*pi*rr*V*d;
%         else
%             f=V*d;
%         end
%         ext_forces(nod(j)*sp-1)=f(1);
%         ext_forces(nod(j)*sp)=f(2);
%     end
% end


% function [Load_nds]=localization(R,T,mats)
% 
%     global GEOMETRY
%     
%     x_0=GEOMETRY.x_0;
%     [nodes,~]=size(x_0);
%     
%     Load_nds=zeros(nodes,mats);
%     for m=1:mats
%         if T(m)==1 || T(m)==4
%             x1=R(1,m*2-1);
%             y1=R(2,m*2-1);
%             d=zeros(nodes,1);
%             for i=1:nodes
%             	d(i)=sqrt((x_0(i,1)-x1)^2+(x_0(i,2)-y1)^2);
%             end
%             [~,j]=min(d);
%             Load_nds(j,m)=1;   
%         elseif T(m)==2 || T(m)==3 || T(m)==5
%             for nodo=1:nodes
%                 tol=GEOMETRY.h_nds(nodo,1)/4;
%                 if (x_0(nodo,2)>=(R(2,m*2-1)-tol)) && (x_0(nodo,2)<=(R(2,m*2)+tol)) 
%                     if (x_0(nodo,1)>=(R(1,m*2-1)-tol)) && (x_0(nodo,1)<=(R(1,m*2)+tol))
%                         Load_nds(nodo,m)=1;
%                     end
%                 end
%             end
%         end
%     end
% 
% end

