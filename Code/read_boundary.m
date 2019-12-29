
function read_boundary(filetxt,BLCK,NODE_LIST)

% File: read_boundary
%   Read loads from boundary.txt
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
    TIED = strings(loads,1);
    INTERVAL=zeros(2,loads); % Intervals of loads
    OUT=zeros(loads,1);
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
                    case 'TIED_NODES'
                        TYPE(M)=5;
                        continue
                    case 'PORE_PRESSURE'
                        TYPE(M)=6;
                        continue
                    otherwise
                        disp('Error, type of boundary not implemented yet!')
                        stop
                end
%             case 'X_RANGE'
%                 cc = str2double(bb{t});
%                 if isnan(cc)
%                     if strcmp(bb{t},'FULL')
%                         RANGE(1,M*2-1)=L;
%                     elseif strcmp(bb{t},'INI')
%                         RANGE(1,M*2-1)=L0;
%                     else
%                         disp('Error, wrong load X range!')
%                     end
%                 else
%                     RANGE(1,M*2-1)=cc;
%                 end
%                 cc = str2double(c{t});
%                 if isnan(cc)
%                     len=strlength(c{t});
%                     if len==0 
%                         RANGE(1,M*2)=RANGE(1,M*2-1);
%                     elseif strcmp(c{t},'FULL')
%                         RANGE(1,M*2)=L;
%                     elseif strcmp(c{t},'INI')
%                         RANGE(1,M*2)=L0;
%                     else
%                         disp('Error, wrong load X range!')
%                         stop
%                     end
%                 else
%                     RANGE(1,M*2)=cc;
%                 end
%                 continue
%             case 'Y_RANGE'
%                 cc = str2double(bb{t});
%                 if isnan(cc)
%                     if strcmp(bb{t},'FULL')
%                         RANGE(2,M*2-1)=H;
%                     elseif strcmp(bb{t},'INI')
%                         RANGE(2,M*2-1)=H0;
%                     else
%                         disp('Error, wrong load Y range!')
%                     end
%                 else
%                     RANGE(2,M*2-1)=cc;
%                 end
%                 cc = str2double(c{t});
%                 if isnan(cc)
%                     len=strlength(c{t});
%                     if len==0 
%                         RANGE(2,M*2)=RANGE(2,M*2-1);
%                     elseif strcmp(c{t},'FULL')
%                         RANGE(2,M*2)=H;
%                     elseif strcmp(c{t},'INI')
%                         RANGE(2,M*2)=H0;
%                     else
%                         disp('Error, wrong load Y range!')
%                         stop
%                     end
%                 else
%                     RANGE(2,M*2)=cc;
%                 end
%                 continue
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
            case 'OUTPUT'
                OUT(M)=str2double(bb{t});
                continue
            case 'NODE_LIST'
                NLIST(M)=bb{t};
                continue
            case 'TIED'
                TIED(M)=bb{t};
                continue
            case 'INTERVAL'
                val=str2double(bb{t});
                val2=str2double(c{t});
                if isnan(val)
                    if strcmp(bb{t},'FULL')
                        if BLCK==1
                            INTERVAL(1,M)=0;
                        else
                            INTERVAL(1,M)=SOLVER.Time_final(BLCK-1);
                        end
                        INTERVAL(2,M)=SOLVER.Time_final(BLCK);
                    elseif strcmp(bb{t},'INI')
                        if BLCK==1
                            INTERVAL(1,M)=0;
                        else
                            INTERVAL(1,M)=SOLVER.Time_final(BLCK-1);
                        end
                        if isnan(val2)
                            if strcmp(c{t},'FULL')
                                INTERVAL(2,M)=SOLVER.Time_final(BLCK);
                            elseif strcmp(c{t},'INI') || strcmp(c{t},'')
                                if BLCK==1
                                    INTERVAL(2,M)=0;
                                else
                                    INTERVAL(2,M)=SOLVER.Time_final(BLCK-1);
                                end
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
            otherwise
            	fprintf('Error, unrecognized sequence: %s !! \n',s1)
                stop
        end
                
    end
    
    fclose(fid);
    
    interval(INTERVAL,VALUE,loads,BLCK);
    %[b_nds]=localization(RANGE,loads,TIED,TYPE);
    
    calculate_boundaries(NLIST,NODE_LIST,VECTOR,TIED,TYPE,loads);
    
    for i=1:loads
        if OUT(i)==1
            tt=0;
            for j=1:length(SOLVER.OutputType(:,1))
                if SOLVER.OutputType(j,1)==0
                    SOLVER.OutputType(j,1)=1;
                    SOLVER.OutputType(j,2)=i;
                    tt=0;
                    break;
                elseif SOLVER.OutputType(j,1)==1 && SOLVER.OutputType(j,2)==i
                    tt=0;
                    break;
                else
                    tt=1;
                end
            end
            if tt==1
                j=j+1;
                SOLVER.OutputType(j,1)=1;
                SOLVER.OutputType(j,2)=i;
            end
        end
    end
end

function interval(INTERVAL,VALUE,loads,BLCK)

    global SOLVER TIME BOUNDARY
    
    if BLCK==1
        ini=1;
        fin=SOLVER.step_final(BLCK);
        BOUNDARY.b_mult=strings(fin,loads);
    else
        ini=SOLVER.step_final(BLCK-1);
        fin=SOLVER.step_final(BLCK);
    end
    
    for m=1:loads
        for i=ini:fin
            t=TIME{BLCK}.t(i);
            if t>=INTERVAL(1,m) && t<=INTERVAL(2,m)
                %BOUNDARY.b_mult(i,m)=1;
                val=str2double(VALUE(m));
                if isnan(val)
                    t=TIME{BLCK}.t(i);
                    val=eval(VALUE(m));
                end
                BOUNDARY.b_mult(i,m)=num2str(val);
            else
                BOUNDARY.b_mult(i,m)='NULL';
            end
        end 
    end

end

function calculate_boundaries(NLIST,NODE_LIST,VECTOR,TIED,TYPE,loads)
    
    global GEOMETRY SOLVER BOUNDARY
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    x=GEOMETRY.x_0;
    
    BOUNDARY.constrains = zeros(GEOMETRY.nodes*df,loads);
    BOUNDARY.dad  = zeros(GEOMETRY.nodes*df,loads);
    BOUNDARY.vad  = zeros(GEOMETRY.nodes*df,loads);
    BOUNDARY.Type = TYPE;
    BOUNDARY.tied = zeros(GEOMETRY.nodes*df,loads); 
    
    BOUNDARY.size=loads;
    
    bcs=NODE_LIST.bcs;
    BC=NODE_LIST.BC;

    for m=1:loads
        % Direction
        V=VECTOR(:,m)';
        nv=norm(V);
        if nv==0
            if TYPE~=6
                continue
            else
                disp('Error on the vector of boundary conditon');
                stop;
            end
        else
            V=V/norm(V);
        end
        
        %NODE LIST
        
        ll=str2double(NLIST(m));
        if isnan(ll)
            if strcmp(NLIST(m),'FULL')
                nod_f=linspace(1,nodes);
            else
                fprintf('Error, unrecognized list of nodes: %s !! \n',ll)
                stop
            end
        else
            if ll>bcs
                fprintf('Error, unrecognized list of nodes: %i !! \n',ll)
                stop
            else
                nod_f=BC{ll};
            end 
        end
        
        if TYPE(m)==5
            l2=str2double(TIED(m));
            if isnan(ll)
                fprintf('Error, unrecognized list of nodes: %s !! \n',ll)
                stop
            end
            if ll>bcs
                fprintf('Error, unrecognized list of nodes: %i !! \n',ll)
                stop
            else
                nod_aux=BC{l2};
            end 
            nod_f2=zeros(length(nod_aux),1);
            for i=1:length(nod_f)
                for j=1:length(nod_aux)
                    if x(nod_f(i),1)==x(nod_aux(j),1) ||...
                            x(nod_f(i),2)==x(nod_aux(j),2) 
                        nod_f2(i)=nod_aux(j);
                        break;
                    end
                end
            end
        end
        
        
        % Nodes
            
        if SOLVER.UW==0 && (TYPE(m)==2 || TYPE(m)==4)
            disp('error, take care of the water boundary conditions!!');
            stop
        elseif SOLVER.UW==1 && TYPE(m)==6
            disp('error, take care of the water boundary conditions!!');
            stop
        elseif SOLVER.UW==2 && (TYPE(m)==2 || TYPE(m)==4)
            disp('error, take care of the water boundary conditions!!');
            stop
        else
            [long,~]=size(nod_f);
            for i=1:long
                if (TYPE(m)==2) || (TYPE(m)==1 && SOLVER.UW==0)
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY.dad(nod_f(i)*df+1-k,m)=V(sp+1-k);
                            BOUNDARY.constrains(nod_f(i)*df+1-k,m)=1;
                        end
                    end
                elseif TYPE(m)==1 && SOLVER.UW==1
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY.dad(nod_f(i)*df-sp+1-k,m)=V(sp+1-k);
                            BOUNDARY.constrains(nod_f(i)*df-sp+1-k,m)=1;
                        end
                    end
                elseif TYPE(m)==1 && SOLVER.UW==2
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY.dad(nod_f(i)*df-k,m)=V(sp+1-k);
                            BOUNDARY.constrains(nod_f(i)*df-k,m)=1;
                        end
                    end
                elseif (TYPE(m)==4) || (TYPE(m)==3 && SOLVER.UW==0)
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY.vad(nod_f(i)*df+1-k,m)=V(sp+1-k);
                            BOUNDARY.constrains(nod_f(i)*df+1-k,m)=2;
                        end
                    end
                elseif TYPE(m)==3 && SOLVER.UW==1
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY.vad(nod_f(i)*df-sp+1-k,m)=V(sp+1-k);
                            BOUNDARY.constrains(nod_f(i)*df-sp+1-k,m)=2;
                        end
                    end
                elseif TYPE(m)==3 && SOLVER.UW==2
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY.vad(nod_f(i)*df-k,m)=V(sp+1-k);
                            BOUNDARY.constrains(nod_f(i)*df-k,m)=2;
                        end
                    end
                elseif TYPE(m)==5
                    for k=1:sp
                        if V(sp+1-k)~=0
                            if SOLVER.UW==1
                                BOUNDARY.dad(nod_f(i)*df-sp+1-k,m)=sign(V(sp+1-k));
                                BOUNDARY.constrains(nod_f(i)*df-sp+1-k,m)=3;
                                BOUNDARY.tied(nod_f(i)*df-sp+1-k,m)=...
                                    nod_f2(i)*df-sp+1-k;
                            elseif SOLVER.UW==2
                                BOUNDARY.dad(nod_f(i)*df-k,m)=sign(V(sp+1-k));
                                BOUNDARY.constrains(nod_f(i)*df-k,m)=3;
                                BOUNDARY.tied(nod_f(i)*df-k,m)=...
                                    nod_f2(i)*df-k;
                            else
                                BOUNDARY.dad(nod_f(i)*df+1-k,m)=sign(V(sp+1-k));
                                BOUNDARY.constrains(nod_f(i)*df+1-k,m)=3;
                                BOUNDARY.tied(nod_f(i)*df+1-k,m)=...
                                    nod_f2(i)*df+1-k;
                            end
                        end
                    end
                elseif TYPE(m)==6 && SOLVER.UW==2
                    BOUNDARY.dad(nod_f(i)*df,m)=1;
                    BOUNDARY.constrains(nod_f(i)*df,m)=1;
                end
                if SOLVER.UW==0 && sp~=df
                    % No water
                    BOUNDARY.constrains(nod_f(i)*df,m)=1;
                    BOUNDARY.constrains(nod_f(i)*df-1,m)=1;
                end
            end
        end    
    end
end

% function [Load_nds]=localization(R,mats,TIED,TYPE)
% 
%     global GEOMETRY
%     x_0=GEOMETRY.x_0;
%     [nodes,~]=size(x_0);
%     Load_nds=zeros(nodes,mats);
%     
%     for m=1:mats
%         for nodo=1:nodes
%             tol=GEOMETRY.h_nds(nodo,1)/5;
%             if (x_0(nodo,2)>=R(2,m*2-1)-tol) && (x_0(nodo,2)<=R(2,m*2)+tol) 
%                 if (x_0(nodo,1)>=R(1,m*2-1)-tol) && (x_0(nodo,1)<=R(1,m*2)+tol)
%                     if TYPE(m)~=5
%                         Load_nds(nodo,m)=1;
%                     else
%                         for k=1:mats
%                             if TYPE(k)==5 && k~=m
%                                 if strcmp(TIED(m),'X')
%                                     D=R(1,k*2-1);
%                                 elseif strcmp(TIED(m),'Y')
%                                     D=R(2,k*2-1);
%                                 else
%                                     disp('no defined tied direction')
%                                     stop
%                                 end 
%                             end
%                         end
%                         for j=1:nodes
%                             if strcmp(TIED(m),'X') && (nodo~=j)
%                                 if (x_0(j,2)>=x_0(nodo,2)-tol) && (x_0(j,2)<=x_0(nodo,2)+tol)
%                                     if (x_0(j,1)>=D-tol) && (x_0(j,1)<=D+tol)
%                                         Load_nds(nodo,m)=j;
%                                         break;
%                                     end
%                                 end
%                             elseif strcmp(TIED(m),'Y') && (nodo~=j)
%                                 if (x_0(j,1)>=x_0(nodo,1)-tol) && (x_0(j,1)<=x_0(nodo,1)+tol)
%                                     if (x_0(j,2)>=D-tol) && (x_0(j,2)<=D+tol) 
%                                         Load_nds(nodo,m)=j;
%                                         break;
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% 
% end
