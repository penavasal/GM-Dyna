
function read_contact

% File: read_contact
%   Read some important problem parameters from contact.txt
%
% Date:
%   Version 1.0   10.07.2018

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some possible parameters if they are not read
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear CONTACT P_SURFACE BODIES
    global CONTACT GEOMETRY P_SURFACE BODIES
    
    x_0=GEOMETRY.x_0;
    [nodes,sp]=size(x_0);
    
    %GEOM
    L=max(x_0(:,1));
    H=max(x_0(:,2));
    
    L0=min(x_0(:,1));
    H0=min(x_0(:,2));
   
    %FILE
    fid = fopen('contact.txt', 'rt'); % opción rt para abrir en modo texto
    formato = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s'; % formato de cada línea 
    data = textscan(fid, formato, 'HeaderLines', 1);
    
    a = data{1};
    % Convertir a vector numérico
    [long,~]=size(a);
    bb = data{2};
    c = data{3};
    d = data{4};
    
    s1=a{1};
    if strcmp(s1,'BODIES')
        bodies=str2double(bb{1});
    else 
        bodies=1;
    end
    
    l=2;
    exit_=0;
    while exit_==0
        s1=a{l};
        if strcmp(s1,'BODY')
            exit_=1;
        elseif strcmp(s1,'CONTACTS') || strcmp(s1,'SURFACES')
            exit_=2;
        end
        l=l+1;
    end
        
    if exit_==2 && bodies==1
        BODIES.node(nodes,1)=0;
        BODIES.element(GEOMETRY.elements,1)=0;
        for nd=1:nodes
            BODIES.node(nd,1)=1;
        end
        for e=1:GEOMETRY.elements
            BODIES.element(e,1)=1;
        end
    elseif exit_==1
        l=l-1;
        while exit_==1
            s1=a{l};
            if strcmp(s1, '//')
                continue
            end
            if strcmp(s1,'BODY')
                if body<=bodies
                    body=str2double(bb{l});
                    continue
                else
                    exit_=2;
                end
            end
            if strcmp(s1,'X_MAX')
                x_max=[str2double(bb{l}) str2double(c{l})];
                continue
            end
            if strcmp(s1,'X_MIN')
                x_min=[str2double(bb{l}) str2double(c{l})];
                for nd=1:nodes
                    tol_b=GEOMETRY.h_nds(nd,1)/4;
                    if (x_0(nd,2)>=(x_min(2)-tol_b)) && (x_0(nd,2)<=(x_max(2)+tol_b)) 
                    if (x_0(nd,1)>=(x_min(1)-tol_b)) && (x_0(nd,1)<=(x_max(1)+tol_b))
                        BODIES.node(nd,1)=body;
                    end
                    end
                end
                for e=1:GEOMETRY.elements
                    tol_b=GEOMETRY.h_ini(e,1)/4;
                    if (GEOMETRY.xg0(e,2)>=(x_min(2)-tol_b)) && ...
                            (GEOMETRY.xg0(e,2)<=(x_max(2)+tol_b)) 
                    if (GEOMETRY.xg0(e,1)>=(x_min(1)-tol_b)) &&...
                            (GEOMETRY.xg0(e,1)<=(x_max(1)+tol_b))
                        BODIES.element(e,1)=body;
                    end
                    end
                end  
                continue
            end
            if strcmp(s1,'CONTACTS') || strcmp(s1,'SURFACES')
                exit_=2;
            end
        end
    end
    
    exit_=0;
    t=l-1;
    while exit_==0
        s1=a{t};    
        if strcmp(s1,'CONTACTS')
            %Allocate
            cons= str2double(bb{t});
            x0  = zeros(cons,sp);
            xc  = zeros(cons,sp);
            tol = zeros(cons,1);
            mu_c= zeros(cons,1);
            MASTER= zeros(cons,1);
            cc  = ones(cons,1);
            CC_ = zeros(10,3);
            for i=1:cons
                CC{i}=CC_;
            end
            % READ
            M=0;
            while exit_==0 && (t<long)       
                t=t+1;
                s1=a{t};
                if strcmp(s1, '//')
                        continue
                end

                if strcmp(s1,'CONTACT')
                    M=str2double(bb{t});
                    if(M>=cons+1)
                        break
                    end
                elseif strcmp(s1,'TOLERANCE')
                    tol(M,1)= str2double(bb{t});
                    continue
                elseif strcmp(s1,'FRICTION')
                    mu_c(M,1)= str2double(bb{t});
                    continue
                elseif strcmp(s1,'MASTER_body')
                    MASTER(M,1)= str2double(bb{t});
                    continue
                elseif strcmp(s1,'X_0')
                    x0(M,1) = str2double(bb{t});
                    x0(M,2) = str2double(c{t});
                    continue
                elseif strcmp(s1,'XC')
                    xc(M,1) = str2double(bb{t});
                    xc(M,2) = str2double(c{t});
                    continue
                elseif strcmp(s1,'CC')
                    CC{M}(cc(M),1)=str2double(bb{t});
                    CC{M}(cc(M),2)=str2double(c{t});
                    CC{M}(cc(M),3)=str2double(d{t});
                    cc(M)=cc(M)+1;
                    continue
                elseif strcmp(s1,'SURFACES')
                    t=t-1;
                    break;
                else    
                    fprintf('Error, unrecognized sequence: %s !! \n',s1)
                    stop
                end   
            end
            if cons
                for m=1:cons
                    if cc(m)<10
                        CC{m}(cc(m):10,:)=[];
                    end
                end
            end
        elseif strcmp(s1,'SURFACES')
            %Allocate
            surfs= str2double(bb{t});
            K    = zeros(surfs,1);
            mu_s = zeros(surfs,1);
            X_   = zeros(10,2);
            xx   = ones(surfs,1);
            for i=1:surfs
                X{i}=X_;
            end
            % READ
            M=0;
            while exit_==0 && (t<long)       
                t=t+1;
                s1=a{t};
                if strcmp(s1, '//')
                        continue
                end

                if strcmp(s1,'POTENTIAL_SURFACE')
                    M=str2double(bb{t});
                    if(M>=surfs+1)
                        break
                    end
                elseif strcmp(s1,'STIFFNESS')
                    K(M,1)= str2double(bb{t});
                    continue
                elseif strcmp(s1,'FRICTION')
                    mu_s(M,1)= str2double(bb{t});
                    continue
                elseif strcmp(s1,'X')
                    X{M}(xx(M),1)=str2double(bb{t});
                    X{M}(xx(M),2)=str2double(c{t});
                    xx(M)=xx(M)+1;
                    continue
                elseif strcmp(s1,'CONTACTS')
                    t=t-1;
                    break;
                else    
                    fprintf('Error, unrecognized sequence: %s !! \n',s1)
                    stop
                end  
            end
            if surfs
                for m=1:surfs
                    if xx(m)<10
                        X{m}(xx(m):10,:)=[];
                    end
                end
            end
        else
            if t<long
                t=t+1;
            else
                exit_=1;
            end
        end
    end
    
    fclose(fid); 
    
    %CONTACT
    if cons
        CONTACT.MASTER=MASTER;
        for  M=1:cons
            [CONTACT.M_LIST{M},CONTACT.nCC(M)]=...
                M_location(x_0,x0(M,:),xc(M,:),CC{M},tol(M,1),MASTER);
        end
    end

    %SURFACE
    if surfs
        P_SURFACE.K=K;
        P_SURFACE.X=X;
        for  M=1:surfs
            ln = length(X{M})-1;
            for i=1:ln
                list(1)=i;
                list(2)=i+1;
                P_SURFACE.PS_LIST{M}(i,:)=list;
                clear list
            end
        end
    end
    
    surfs;
end

function [LIST,k]=M_location(x_a,x0,xc,CC,tol,MASTER)

    global GEOMETRY BODIES

    [nCC,~]=size(CC);

    k=0;
    for i=1:nCC
        l=0;
        clear list
        if CC(i,1)==4
            R=CC(i,2);
            A=CC(i,3);
            [list2]=circ(xc,x0,x_a,R,A,tol);
            for j=1:length(list2)-1
                ll=[list2(j); list2(j+1)];
                k=k+1;
                LIST{k}=ll;
            end 
            x0=x_a(list2(length(list2)),:);
        else
            if CC(i,1)==1
                for j=1:GEOMETRY.nodes
                    d=x0(1)-x_a(j,1);
                    if abs(d)<tol && BODIES.node(j)==MASTER
                        l=l+1;
                        list(l)=j;
                    end               
                end 
            else

                M=CC(i,2);
                rM=sqrt(M^2+1);
                for j=1:GEOMETRY.nodes
                    d=(M*x_a(j,1)-x_a(j,2)+x0(1,2)-M*x0(1,1))/rM;
                    if abs(d)<tol && BODIES.node(j)==MASTER
                        l=l+1;
                        list(l)=j;
                    end               
                end       
            end 
            [list]=order(list,x0,x_a);
            x0=x_a(list(length(list)),:);
            k=k+1;
            LIST{k}=list;
        end
    end

end


function [list2]=circ(xc,x0,x_a,R,A,tol)

    [nodes,~]=size(x_a);
    
    %%%Select
    l=0;
    for i=1:nodes
        d=sqrt((x_a(i,1)-xc(1))^2+(x_a(i,2)-xc(2))^2);
        if abs(d)<R+tol  && abs(d)>R-tol
            l=l+1;
            list(l)=i;
        end 
    end
    
    %%%Order
    
    %%% 1st and 2nd
    for i=1:l
        if x0(1)==x_a(list(i),1) && x0(2)==x_a(list(i),2)
            list2(1)=list(i);
            list3(1)=i;
            break;
        end
    end
    [j,k]=dist2(x_a,list,i);
    
    V=[x_a(list2(1),1)-xc(1) x_a(list2(1),2)-xc(2)];
    
    P_vec=(x_a(list2(1),1)-xc(1))*(x_a(list(j),2)-xc(2))-...
        (x_a(list(j),1)-xc(1))*(x_a(list2(1),2)-xc(2));
    
    if P_vec>0
        list2(2)=list(j);
        list3(2)=j;
        V1=[x_a(list(j),1)-xc(1) x_a(list(j),2)-xc(2)];
    else
        list2(2)=list(k);
        list3(2)=k;
        V1=[x_a(list(k),1)-xc(1) x_a(list(k),2)-xc(2)];
    end
    
    i=2;
    Angle_old=acos((V(1)*V1(1)+V(2)*V1(2))/...
            sqrt(V(1)^2+V(2)^2)/sqrt(V1(1)^2+V1(2)^2));
    while i
        [j]=dist(x_a,list,list3(i),list3(i-1));
        V1=[x_a(list(j),1)-xc(1) x_a(list(j),2)-xc(2)];
        Angle=acos((V(1)*V1(1)+V(2)*V1(2))/...
            sqrt(V(1)^2+V(2)^2)/sqrt(V1(1)^2+V1(2)^2));
        %Exact angle
        if Angle<Angle_old
            if Angle_old>=pi
                Angle=2*pi-Angle;
            end
        end
        
        %Final point
        if Angle<A+tol && Angle>A-tol
            i=i+1;
            list2(i)=list(j);
            list3(i)=j;
            i=0;
        elseif Angle<A
            i=i+1;
            list2(i)=list(j);
            list3(i)=j;
            Angle_old=Angle;
        else
            i=0;
        end
    end    
end

function [i]=dist(x_a,nd,j,k)
    x1=x_a(nd(j),1);
    y1=x_a(nd(j),2);
    d=zeros(length(nd),1);
    for i=1:length(nd)
        if(i~=j) && (i~=k)
            d(i)=sqrt((x_a(nd(i),1)-x1)^2+(x_a(nd(i),2)-y1)^2);
        else
            d(i)=1.0e32;
        end
    end
    [~,i]=min(d);  
end

function [i,k]=dist2(x_a,nd,j)
    x1=x_a(nd(j),1);
    y1=x_a(nd(j),2);
    d=zeros(length(nd),1);
    for i=1:length(nd)
        if(i~=j)
            d(i)=sqrt((x_a(nd(i),1)-x1)^2+(x_a(nd(i),2)-y1)^2);
        else
            d(i)=1.0e32;
        end
    end
    [~,i]=min(d);
    d(i)=1.0e32;
    [~,k]=min(d);  
end

function [l2]=order(l1,x0,x_a)

    N=length(l1);

    d=zeros(N,1);
    for i=1:N
        d(i)=sqrt((x0(1,1)-x_a(l1(i),1))^2+(x0(1,2)-x_a(l1(i),2))^2);
    end

    l2=zeros(N,1);
    for i=1:N
        [~,b]=min(d);
        l2(i)=l1(b);
        d(b)=1e32;
    end
    N;
end
