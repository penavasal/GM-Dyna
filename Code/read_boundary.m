
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

    fid = fopen(filetxt, 'rt'); % opci�n rt para abrir en modo texto
    formato = '%s %s %s %s'; % formato de cada l�nea 
    data = textscan(fid, formato, 'HeaderLines', 1);
    
    a = data{1};
    % Convertir a vector num�rico
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
    DIST = strings(loads,1);
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
                    case 'ABSORBING_BC'
                        TYPE(M)=7;
                        SOLVER.absorbing = 1;
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
            case 'DISTRIBUTION'
                DIST(M)=bb{t};
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
    
    % BC for Mixed element
    if (strcmp(GEOMETRY.ELEMENT,'Q8P4') || strcmp(GEOMETRY.ELEMENT,'Q8P4-4'))...
        || (strcmp(GEOMETRY.ELEMENT,'T6P3') || strcmp(GEOMETRY.ELEMENT,'T6P3-3'))
        loads=loads+1;
        NLIST(loads)='Q8';
        VALUE(loads)='0';
        if BLCK>1
            INTERVAL(1,loads)=SOLVER.Time_final(BLCK-1);
        else
            INTERVAL(1,loads)=0;
        end
        INTERVAL(2,loads)=SOLVER.Time_final(BLCK);
        TIED(loads)=0;
        OUT(loads)=0;
        DIST(loads)="";
        if SOLVER.UW==1 || SOLVER.UW==4
            TYPE(loads)=2;
            VECTOR(:,loads)=[1; 1];
        elseif SOLVER.UW==2 || SOLVER.UW==3
            TYPE(loads)=6;
            VECTOR(:,loads)=[0; 0];
        end
    end
    
    
    interval(INTERVAL,VALUE,loads,BLCK);
    %[b_nds]=localization(RANGE,loads,TIED,TYPE);
    
    calculate_boundaries(NLIST,NODE_LIST,VECTOR,TIED,TYPE,loads,BLCK,DIST);
   
    
    % OUTPUT
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

    global BOUNDARY
    
    BOUNDARY{BLCK}.b_mult=strings(4,loads);
    %BOUNDARY{BLCK}.b_mult(4,loads)=0;
    for m=1:loads
        % Value
        val=str2double(VALUE(m));
        if isnan(val)
            BOUNDARY{BLCK}.b_mult(4,m)='FUNCTION';
            BOUNDARY{BLCK}.b_mult(1,m)=VALUE(m);
        else
            BOUNDARY{BLCK}.b_mult(4,m)='VALUE';
            BOUNDARY{BLCK}.b_mult(1,m)=val;
        end
        BOUNDARY{BLCK}.b_mult(2,m)=INTERVAL(1,m);
        BOUNDARY{BLCK}.b_mult(3,m)=INTERVAL(2,m);
    end

end

function calculate_boundaries(NLIST,NODE_LIST,VECTOR,TIED,TYPE,loads,BLCK,DIST)
    
    global GEOMETRY SOLVER BOUNDARY
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    x0=GEOMETRY.x_0;
    
    BOUNDARY{BLCK}.constrains = zeros(GEOMETRY.nodes*df,loads);
    BOUNDARY{BLCK}.dad  = zeros(GEOMETRY.nodes*df,loads);
    BOUNDARY{BLCK}.vad  = zeros(GEOMETRY.nodes*df,loads);
    BOUNDARY{BLCK}.Type = TYPE;
    BOUNDARY{BLCK}.tied = zeros(GEOMETRY.nodes*df,loads); 
    BOUNDARY{BLCK}.abc = {};
    
    BOUNDARY{BLCK}.size=loads;
    
    bcs=NODE_LIST.bcs;
    BC=NODE_LIST.BC;
    abcs=NODE_LIST.abcs;
    ABC=NODE_LIST.ABC;
    
   
    for m=1:loads
        
        dd=DIST(m);
        rr=strlength(dd);
        if rr==0
            dd='1';
        end
        
        
        % Direction
        V=VECTOR(:,m)';
        nv=norm(V);
        if nv==0
            if TYPE(m)~=6 && TYPE(m)~=7
                disp('Error on the vector of boundary conditon');
                stop;
            end
        else
            V=V/norm(V);
        end
        
        %NODE LIST
        
        if TYPE(m)<4 || TYPE(m)==6
            ll=str2double(NLIST(m));
            if isnan(ll)
                if strcmp(NLIST(m),'FULL')
                    nod_f=linspace(1,GEOMETRY.nodes,GEOMETRY.nodes)';
                elseif strcmp(NLIST(m),'Q8')
                    nod_f=setdiff(GEOMETRY.elem,GEOMETRY.elem_c,'sorted');
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
        elseif TYPE(m)==7
            ll=str2double(NLIST(m));
            if isnan(ll)
                fprintf('Error, unrecognized list of nodes: %s !! \n',ll)
                stop
            else
                if ll>abcs
                    fprintf('Error, unrecognized list of nodes: %i !! \n',ll)
                    stop
                else
                    nod_f=ABC{ll};
                end 
            end
            
            BOUNDARY{BLCK}.abc{1,m}=nod_f;
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
                    if x0(nod_f(i),1)==x0(nod_aux(j),1) ||...
                            x0(nod_f(i),2)==x0(nod_aux(j),2) 
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
        elseif (SOLVER.UW==1 || SOLVER.UW==4)&& TYPE(m)==6
            disp('error, take care of the water boundary conditions!!');
            stop
        elseif SOLVER.UW==2 && (TYPE(m)==2 || TYPE(m)==4)
            disp('error, take care of the water boundary conditions!!');
            stop
        elseif TYPE(m)~=7
            [long,~]=size(nod_f);
            for i=1:long
                
                nn=nod_f(i);
                if ~strcmp(dd,'1')
                    x=x0(nn,1);
                    y=x0(nn,2);
                end
                d=eval(dd);
                
                
                if (TYPE(m)==2 && SOLVER.UW==1) || (TYPE(m)==1 && SOLVER.UW==0)
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY{BLCK}.dad(nod_f(i)*df+1-k,m)=d*V(sp+1-k);
                            BOUNDARY{BLCK}.constrains(nod_f(i)*df+1-k,m)=1;
                        end
                    end
                elseif TYPE(m)==1 && (SOLVER.UW==1 || SOLVER.UW==4)
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY{BLCK}.dad(nod_f(i)*df-sp+1-k,m)=d*V(sp+1-k);
                            BOUNDARY{BLCK}.constrains(nod_f(i)*df-sp+1-k,m)=1;
                        end
                    end
                elseif (TYPE(m)==1 && SOLVER.UW==2) || (TYPE(m)==2 && SOLVER.UW==3)
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY{BLCK}.dad(nod_f(i)*df-k,m)=d*V(sp+1-k);
                            BOUNDARY{BLCK}.constrains(nod_f(i)*df-k,m)=1;
                        end
                    end
                elseif TYPE(m)==1 && SOLVER.UW==3
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY{BLCK}.dad(nod_f(i)*df-sp-k,m)=d*V(sp+1-k);
                            BOUNDARY{BLCK}.constrains(nod_f(i)*df-sp-k,m)=1;
                        end
                    end
                elseif (TYPE(m)==4) || (TYPE(m)==3 && SOLVER.UW==0)
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY{BLCK}.vad(nod_f(i)*df+1-k,m)=d*V(sp+1-k);
                            BOUNDARY{BLCK}.constrains(nod_f(i)*df+1-k,m)=2;
                        end
                    end
                elseif TYPE(m)==3 && (SOLVER.UW==1 || SOLVER.UW==4)
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY{BLCK}.vad(nod_f(i)*df-sp+1-k,m)=d*V(sp+1-k);
                            BOUNDARY{BLCK}.constrains(nod_f(i)*df-sp+1-k,m)=2;
                        end
                    end
                elseif (TYPE(m)==3 && SOLVER.UW==2) || (TYPE(m)==4 && SOLVER.UW==3)
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY{BLCK}.vad(nod_f(i)*df-k,m)=d*V(sp+1-k);
                            BOUNDARY{BLCK}.constrains(nod_f(i)*df-k,m)=2;
                        end
                    end
                elseif TYPE(m)==3 && SOLVER.UW==3 
                    for k=1:sp
                        if V(sp+1-k)~=0
                            BOUNDARY{BLCK}.vad(nod_f(i)*df-sp-k,m)=d*V(sp+1-k);
                            BOUNDARY{BLCK}.constrains(nod_f(i)*df-sp-k,m)=2;
                        end
                    end
                elseif TYPE(m)==5
                    for k=1:sp
                        if V(sp+1-k)~=0
                            if SOLVER.UW==1 || SOLVER.UW==4
                                BOUNDARY{BLCK}.dad(nod_f(i)*df-sp+1-k,m)=d*sign(V(sp+1-k));
                                BOUNDARY{BLCK}.constrains(nod_f(i)*df-sp+1-k,m)=3;
                                BOUNDARY{BLCK}.tied(nod_f(i)*df-sp+1-k,m)=...
                                    nod_f2(i)*df-sp+1-k;
                            elseif SOLVER.UW==3
                                BOUNDARY{BLCK}.dad(nod_f(i)*df-sp-k,m)=d*sign(V(sp+1-k));
                                BOUNDARY{BLCK}.constrains(nod_f(i)*df-sp-k,m)=3;
                                BOUNDARY{BLCK}.tied(nod_f(i)*df-sp-k,m)=...
                                    nod_f2(i)*df-sp-k;
                            elseif SOLVER.UW==2
                                BOUNDARY{BLCK}.dad(nod_f(i)*df-k,m)=d*sign(V(sp+1-k));
                                BOUNDARY{BLCK}.constrains(nod_f(i)*df-k,m)=3;
                                BOUNDARY{BLCK}.tied(nod_f(i)*df-k,m)=...
                                    nod_f2(i)*df-k;
                            else
                                BOUNDARY{BLCK}.dad(nod_f(i)*df+1-k,m)=d*sign(V(sp+1-k));
                                BOUNDARY{BLCK}.constrains(nod_f(i)*df+1-k,m)=3;
                                BOUNDARY{BLCK}.tied(nod_f(i)*df+1-k,m)=...
                                    nod_f2(i)*df+1-k;
                            end
                        end
                    end
                elseif TYPE(m)==6 && (SOLVER.UW==2 || SOLVER.UW==3)
                    BOUNDARY{BLCK}.dad(nod_f(i)*df,m)=1;
                    BOUNDARY{BLCK}.constrains(nod_f(i)*df,m)=1;
                end
                if SOLVER.UW==0 && sp~=df
                    % No water
                    BOUNDARY{BLCK}.constrains(nod_f(i)*df,m)=1;
                    BOUNDARY{BLCK}.constrains(nod_f(i)*df-1,m)=1;
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
