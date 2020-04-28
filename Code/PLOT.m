function PLOT(driver,folder,varargin)

    ll=length(varargin);
    b = mod(ll+1,2);
    a = floor((ll)/2);
    
    %Initial values
    ampl=1;
    freq=1;
    steps=0;
    str_out='RES';
    PLT=struct('var',0,'film',0,'near',0,'file',0,'e',0,'range',0,...
        'onnodes',0,'time',0,'pw',1,'drained',0);

    if ~strcmp(driver,'NEIGHBORS') && b~=0
        error('Wrong format: file,type,num,type,num...');
    elseif strcmp(driver,'NEIGHBORS') && mod(nargin,2)~=0
        error('Wrong format: file,type,num,type,num...');
    else
        t=0;
        if ~strcmp(driver,'NEIGHBORS')
            t=1;
            str_in=varargin{t};
        end
        
        for i=1:a
            t=t+1;
            s2=varargin{t};
            switch s2
                case 'AMPL'
                    t=t+1;
                    ampl=varargin{t};
                    continue
                case 'AMPLIFICATION'
                    t=t+1;
                    ampl=varargin{t};
                    continue
                case 'FREQ'
                    t=t+1;
                    freq=varargin{t};
                    continue
                case 'FREQUENCY'
                    t=t+1;
                    freq=varargin{t};
                    continue
                case 'FILE_OUT'
                    t=t+1;
                    str_out=varargin{t};
                    continue
                case 'OUT'
                    t=t+1;
                    str_out=varargin{t};
                    continue
                case 'VAR'
                    t=t+1;
                    PLT.var=varargin{t};
                    continue
                case 'FILE_IN'
                    t=t+1;
                    PLT.file=varargin{t};
                    continue
                case 'E'
                    t=t+1;
                    PLT.e=varargin{t};
                    continue
                case 'FILM'
                    t=t+1;
                    PLT.film=1;
                    continue
                case 'TIME'
                    t=t+1;
                    PLT.time=varargin{t};
                    continue
                case 'RANGE'
                    t=t+1;
                    PLT.range=varargin{t};
                    continue
                case 'ON_NODES'
                    t=t+1;
                    PLT.onnodes=1;
                    continue
                case 'PW'
                    t=t+1;
                    PLT.pw=varargin{t};
                    continue
                case 'DRAINED'
                    t=t+1;
                    PLT.drained=1;
                    continue
                case 'STEPS'
                    t=t+1;
                    if strcmp(varargin{t},'FULL')
                        steps=0;
                    else
                        steps=varargin{t};
                    end
                    continue
                otherwise
                    error('Unrecognized variable')
            end
        end

    end

    if strcmp(driver,'VTK')
        vtk_driver(str_in,str_out,steps,ampl,freq,folder)
    elseif strcmp(driver,'GID')
        GID_driver(str_in,str_out,steps,freq,folder)
    elseif strcmp(driver,'GRID')
        plot_ep(str_in,folder,steps,freq,PLT)
    elseif strcmp(driver,'DEFORMED') || strcmp(driver,'DISTRIBUTION')
        plot_f(driver,str_in,folder,steps,ampl,freq,PLT)
    elseif strcmp(driver,'CRACK')
        plot_mp(str_in,folder,steps,freq)
    elseif strcmp(driver,'NEIGHBORS')
        plot_nb(PLT,ampl)
    elseif strcmp(driver,'CONSTITUTIVE')
        plot_constitutive(str_in,steps,freq,PLT)
    elseif strcmp(driver,'EXCESS_PW')
        nodes_col(str_in,PLT)
    elseif strcmp(driver,'CONSOLIDATION')
        plot_col(str_in,PLT)
    end
    
end

function plot_nb(PLT,h)

if PLT.file==0
    global GEOMETRY
    ds=zeros(GEOMETRY.nodes*GEOMETRY.df,1);
else
    load(PLT.file,'-mat','GEOMETRY','GLOBAL');
    ds=GLOBAL.d(:,GLOBAL.ste_p);
end

NODO=PLT.e;
nds=PLT.near;

xg=GEOMETRY.xg_0;
x_a=GEOMETRY.x_0;
elem=GEOMETRY.elem_c;
df=GEOMETRY.df;
nodes=GEOMETRY.nodes;

[~,NNE]=size(elem);

if h~=0
    for i=1:nodes
        x_a(i,1)=x_a(i,1)+h*ds((i-1)*df+1,1);
        x_a(i,2)=x_a(i,2)+h*ds((i-1)*df+2,1);
    end
end


if NODO
    XX=zeros(length(nds),2);
    for i=1:length(nds)
        XX(i,1)=x_a(nds(i),1);
        XX(i,2)=x_a(nds(i),2);
    end
    if NNE==3
    hold on, triplot(elem,x_a(:,1),x_a(:,2)), scatter(xg(:,1),xg(:,2)), ...
        scatter(xg(NODO,1),xg(NODO,2),'MarkerFaceColor',[1 1 0]), ...
        scatter(XX(:,1),XX(:,2),'MarkerFaceColor',[1 0 0]), hold off
    elseif NNE==4
    hold on, quadplot(elem,x_a(:,1),x_a(:,2)), scatter(xg(:,1),xg(:,2)), ...
        scatter(xg(NODO,1),xg(NODO,2),'MarkerFaceColor',[1 1 0]), ...
        scatter(XX(:,1),XX(:,2),'MarkerFaceColor',[1 0 0]), hold off
    end
else
    if NNE==3
    hold on, triplot(elem,x_a(:,1),x_a(:,2)), scatter(xg(:,1),xg(:,2)), ...
        hold off
    elseif NNE==4
    hold on, quadplot(elem,x_a(:,1),x_a(:,2)), scatter(xg(:,1),xg(:,2)), ...
        hold off
    end
    xlabel('X'),ylabel('Y')
end
    

end

function plot_mp(str,folder,steps,rr)

    root=strcat(folder,'/',str,'.mat');
    load(root,'GEOMETRY','SOLVER','GLOBAL');

    Status=GLOBAL.status;
    d = GLOBAL.d;
    x_0=GEOMETRY.x_0;
    xg=GLOBAL.xg;
    Elements=GEOMETRY.mat_points;
    sp=GEOMETRY.sp;
    Nodos=GEOMETRY.nodes;
    
    if steps==0
        contador=GLOBAL.ste_p-1;
    else
        contador=min(steps,ste_p);
    end

    DDD=max(x_0(:,1));
    HHH=max(x_0(:,2));

    %
    for cont=1:rr:contador
        s1=0;
        s=0;
        for nodo=1:Nodos
            if GEOMETRY.body(1)==1
                s=s+1;
                x2(s,1)=x_0(nodo,1)+d(sp*nodo-1,cont);
                y2(s,1)=x_0(nodo,2)+d(sp*nodo,cont);
            else
                s1=s1+1;
                x1(s1,1)=x_0(nodo,1)+d(sp*nodo-1,cont);
                y1(s1,1)=x_0(nodo,2)+d(sp*nodo,cont);
            end
        end
        t=0;
        for e=1:Elements
            if Status(e,cont)==1
                t=t+1;
                x3(t,1)=xg(e,cont);
                y3(t,1)=xg(Elements+e,cont);
            end
        end
        
        scatter(x2,y2,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]), ...
        hold on,...
        if s1
            scatter(x1,y1,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]), ...
        end
        if t
            scatter(x3,y3,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'Marker','.'),...
        end
        hold off
        axis([-DDD/10,DDD+DDD/10,-HHH/10,HHH+HHH/10])
        drawnow
        
        clear x2 y2 y1 x1 x3 y3
    end
    

end

function plot_ep(str,folder,ste_p,rt,PLT)

    
    root=strcat(folder,'/',str,'.mat');
    load(root,'GEOMETRY','SOLVER','GLOBAL');
    
    if ste_p==0
        ste_p=GLOBAL.ste_p-1;
    end
    
    %[elements,NNE]=size(GEOMETRY.elem);
    x_0=GEOMETRY.x_0;
    
    if PLT.var==0
        error('Unrecognized variable')
    else
        if strcmp(PLT.var,'GAMMA')
            VAR=GLOBAL.gamma_nds;
        %elseif strcmp(varargin{1},'PW')
        else
            error('Unrecognized variable')
        end
    end

    [Nodos,~]=size(x_0);
    x=zeros(Nodos,1);
    y=zeros(Nodos,1);
    for nd=1:Nodos
        x(nd,1)=x_0(nd,1);
        y(nd,1)=x_0(nd,2);
    end

    DDD=max(x(:));
    HHH=max(y(:));

    figure
    [xg1,yg1]=meshgrid(0:0.05:DDD,0:0.1:HHH);

    for cont=1:rt:ste_p
        PWW=griddata(x,y,VAR(:,cont),xg1,yg1,'cubic');
        surf(xg1,yg1,PWW)
        colorbar
        axis([0,DDD,0,HHH])
        drawnow
    end

    %Plot contour of pore pressure to obtain 0 pore pressure value
    figure
    %for cont=1:5:ste_p
        PW=griddata(x,y,VAR(:,cont),xg1,yg1,'cubic');
        [c,h]=contour('v6',xg1,yg1,PW);
        clabel(c,h);
        %cont
    %end

end

function plot_f(driver,str,folder,ste_p,hh,rr,PLT)

    %hh   %Amplification in pressure plot
    %rr   %Frequency in pressure plot
    film=0; %film or not?
    film2=0;
    
    flag2=0; %Displacements of the whole
    flag3=0; %Spatial distribution
    
    % READ
    root=strcat(folder,'/',str,'.mat');
    load(root,'GEOMETRY','SOLVER','GLOBAL');


    %TYPE
    if strcmp(driver,'DISTRIBUTION')
        flag3=1;
        if PLT.var==0
            error('Unrecognized variable')
        else
            if strcmp(PLT.var,'PS')
                VAR=GLOBAL.Ps;
            elseif strcmp(PLT.var,'PW')
                VAR=GLOBAL.pw;
            elseif strcmp(PLT.var,'GAMMA')
                VAR=GLOBAL.gamma;
            else
                error('Unrecognized variable')
            end
        end
        if PLT.film
            film2=1;
        end  
    elseif strcmp(driver,'DEFORMED')
        flag2=1;
        if PLT.film
            film=1;
        end
    end

    % Steps
    if ste_p==0
        ste_p=GLOBAL.ste_p-1;
    end
    contador=ste_p;
    
    
    % Allocate
    
    [elements,NNE]=size(GEOMETRY.elem);
     x_0=GEOMETRY.x_0;
     df=GEOMETRY.df;
    
    if NNE==4 || NNE==3
        elem=GEOMETRY.elem;
        elem2=GEOMETRY.elem;
        nodes_p = linspace(1,GEOMETRY.nodes,GEOMETRY.nodes)';
    elseif NNE==8
        NNE=4;
        elem=GEOMETRY.elem_c;
        nodes_p = intersect(GEOMETRY.elem,GEOMETRY.elem_c);
        elem2=zeros(elements,4);
        for i=1:4
            for j=1:elements
                for k=1:length(nodes_p)
                    if elem(j,i)==nodes_p(k)
                        elem2(j,i)=k;
                        break;
                    end
                end
            end
        end
    elseif NNE==6
        NNE=3;
        elem=GEOMETRY.elem_c;
        nodes_p = intersect(GEOMETRY.elem,GEOMETRY.elem_c);
        elem2=zeros(elements,3);
        for i=1:3
            for j=1:elements
                for k=1:length(nodes_p)
                    if elem(j,i)==nodes_p(k)
                        elem2(j,i)=k;
                        break;
                    end
                end
            end
        end
    end
    


    Nodos=length(nodes_p);
    x=zeros(Nodos,1);
    y=zeros(Nodos,1);

    for j=1:Nodos
        nd=nodes_p(j);
        x(j,1)=x_0(nd,1);
        y(j,1)=x_0(nd,2);
    end

    DDD=max(x(:));
    HHH=max(y(:));


    %
    if flag2==1

        if film
            movie = VideoWriter('video.avi');
            open(movie);
        end

        x2=zeros(Nodos,1);
        y2=zeros(Nodos,1);
        figure
        for cont=1:rr:contador
            for j=1:Nodos
                nodo=nodes_p(j);
                if SOLVER.UW==1
                    x2(j)=x(j)+hh*GLOBAL.d(df*nodo-3,cont);
                    y2(j)=y(j)+hh*GLOBAL.d(df*nodo-2,cont);
                elseif SOLVER.UW==2
                    x2(j)=x(j)+hh*GLOBAL.d(df*nodo-2,cont);
                    y2(j)=y(j)+hh*GLOBAL.d(df*nodo-1,cont);
                else
                    x2(j)=x(j)+hh*GLOBAL.d(df*nodo-1,cont);
                    y2(j)=y(j)+hh*GLOBAL.d(df*nodo,cont);
                end
            end
            if NNE==3
                triplot(elem2,x2,y2)
            elseif NNE==4
                quadplot(elem2,x2,y2)
            end
            axis([-0.05,DDD*1.2,-0.05,HHH*1.1])
            drawnow

            if film
                   frame = getframe;
                   writeVideo(movie,frame);  
            end
        end

        if film
            close(movie);
        end
    end

    % Pressure evolution 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag3

        if film2
            movie2 = VideoWriter('press.avi');
            open(movie2);
        end

        PWnodo=zeros(Nodos,contador);
        x1=zeros(Nodos,1);
        y1=zeros(Nodos,1);
        % Nodal Pw2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        for nodo=1:Nodos
            for cont=1:contador
                suma1=0;
                suma2=0;
                for i=1:elements
                    for j=1:NNE
                        if nodo==elem2(i,j)
                            suma1=suma1+VAR(i,cont)*GEOMETRY.Area(i);
                            suma2=suma2+GEOMETRY.Area(i);
                        end
                    end
                end
                PWnodo(nodo,cont)=suma1/suma2;
            end
        end


        figure
        for cont=1:rr:ste_p       
            for j=1:Nodos
                nodo=nodes_p(j);
                if SOLVER.UW==1
                    x1(j)=x(j)+hh*GLOBAL.d(df*nodo-3,cont);
                    y1(j)=y(j)+hh*GLOBAL.d(df*nodo-2,cont);
                elseif SOLVER.UW==2
                    x1(j)=x(j)+hh*GLOBAL.d(df*nodo-2,cont);
                    y1(j)=y(j)+hh*GLOBAL.d(df*nodo-1,cont);
                else
                    x1(j)=x(j)+hh*GLOBAL.d(df*nodo-1,cont);
                    y1(j)=y(j)+hh*GLOBAL.d(df*nodo,cont);
                end
            end
            pr=TriRep(elem2,x1(:),y1(:),PWnodo(:,cont));
            trisurf(pr)
            colorbar
            axis([-0.05,DDD*1.2,-0.05,HHH*1.1])
            drawnow

            if film2
                   frame = getframe;
                   writeVideo(movie2,frame);  
            end
        end

        if film2
            close(movie2);
        end
    end


end

function plot_constitutive(str,steps,freq,PLT)

global GEOMETRY MATERIAL SOLVER

time_d=0.001;

e=PLT.e;
dr=PLT.drained;

load(str,'-mat','GLOBAL','GEOMETRY','MATERIAL','SOLVER');

BLCK=1;

mati=GEOMETRY.material(e);
MODEL=MATERIAL(BLCK).MODEL(mati,1);
MAT=MATERIAL(BLCK).MAT;

if MODEL>=2
    if ~isempty(MAT{19,mati})
        M = MAT{19,mati};
    else
        M=0;
    end
end

if ~isempty(MAT{16,mati})
    n0=MAT{16,mati};
    %n=1-(1-n0)/MAT_POINT{1}(e).J;
else
    n0=0;
end

e_0=n0/(1-n0);

Ps=GLOBAL.Ps(e,:);
P0=Ps(1);

if P0<500 % To kPa
    mult=1000;
else
    mult=0.001;
end
Ps=Ps*mult;
Sy_tot=GLOBAL.Sy(e,:)*mult;
Qs=GLOBAL.Qs(e,:)*mult;
Pw=GLOBAL.pw(e,:)*mult;

Es=GLOBAL.Es;
Es_p=GLOBAL.Es_p;
void_index=GLOBAL.J(e,:)*(1+e_0)-1;


% Time
if steps==0
    ste_p=GLOBAL.ste_p-1;
else
    ste_p=steps;
end

if freq==0
    total=100;
    if total>ste_p
        each=1;
    else
        each=round(ste_p/total);
    end
else
    each=freq;
end
    

% Epsilon
edev=zeros(ste_p,1);
evol=zeros(ste_p,1);
for j=1:ste_p
    [evol(j),edev(j)]=VECTORS.E_invar(Es(:,j)+Es_p(:,j),e);
end

%Ellipse

if MODEL>=3 && MODEL<4
    Pc_max=max(-Sy_tot(e,1:ste_p));
else
    Pc_max=max(-Ps(e,1:ste_p))*1.2;
end
b_max=M*Pc_max/2;
lim=1.2*Pc_max;

figure;

hold on
h2=plot(0,0);
h3=plot(0,0);
for i=1:each:ste_p
    
    if  MODEL>=2 && MODEL<5
        if MODEL>=2 && MODEL<3
            C=-Sy_tot(i);
        else
            C=0;
        end
        delete(h3)
        subplot(2,2,2)
        hold on
        h3=plot(linspace(0,lim,5),M*linspace(0,lim,5)+C,'k');
    else
        C=0;
    end
    
    subplot(2,2,1)
    if dr==0
        plot(edev(1:i)*100,Qs(1:i))
        axis([0 inf 0 max(max(Qs(1:ste_p))*1.1,max(b_max,lim*M+C))])
        xlabel('\epsilon_s %')
        ylabel('Q [kPa]')
    else
        plot(edev(1:i)*100,Qs(e,1:i)./Ps(e,1:i))
        xlabel('\epsilon_s %')
        ylabel('\eta')
    end

    subplot(2,2,2)
    axis([0 Pc_max 0 max(max(Qs(1:ste_p))*1.1,max(b_max,lim*M+C))])
    xlabel('P [kPa]')
    ylabel('Q [kPa]')

    if MODEL>=3 && MODEL<4
    delete(h2)
    Pc=-Sy_tot(i);
    b=M*Pc/2;
    x0=Pc/2;
    t=0:0.01:pi;
    x=x0+Pc/2*cos(t);
    y=b*sin(t);
    
    h2=plot(x,y,'b');
    hold on
    end
    
    plot(-Ps(1:i),Qs(1:i),'r')
    
    if dr==0
        if SOLVER.UW
        subplot(2,2,3)
        plot(edev(1:i)*100,Pw(1:i))
        axis([0 inf -inf max(Pw(1:ste_p))*1.1])
        xlabel('\epsilon %')
        ylabel('P_w [kPa]')
        end
    else
        subplot(2,2,3)
        plot(edev(1:i)*100,evol(1:i)*100)
        xlabel('\epsilon_s %')
        ylabel('\epsilon_v %')
    end

    subplot(2,2,4)
    if dr==0
        plot(-Ps(2:i),1+void_index(2:i))
        axis([0 Pc_max min(void_index(2:ste_p))+0.5 max(void_index(2:ste_p))+1.5])
        xlabel('P [kPa]')
        ylabel('1+e')
    else
        plot(Qs(e,1:i)./Ps(e,1:i),evol(1:i)*100)
        xlabel('\eta')
        ylabel('\epsilon_v %')
    end
    
    drawnow;
    pause(time_d)
end
hold off

end

function nodes_col(str,PLT)

    load(str,'-mat','GLOBAL','GEOMETRY')
    
    rr=PLT.range;
    range1=rr(1);
    range2=rr(2);
    tt=PLT.time;
    
    %nds=GEOMETRY.nodes;
    mps=GEOMETRY.mat_points;
    %sp=GEOMETRY.sp;
    %df=GEOMETRY.df;
    %xa=GEOMETRY.x_0;
    xg=GEOMETRY.xg_0;
    %d=GLOBAL.d;
    
    Pw=GLOBAL.pw;
    S=GLOBAL.Sigma;
    
    tp=GLOBAL.tp;
    steps=GLOBAL.ste_p;
    
    k=0;
    for i=1:mps
        if xg(i,1)<range2  && xg(i,1)>range1
            k=k+1;
            col(k,1)=i;
            col(k,2)=xg(i,2);
        end
    end

    column=zeros(k,2);
    s0=zeros(k,1);
    for i=1:k
        [~,j]=max(col(:,2));
        column(i,1)=col(j,1);
        column(i,2)=xg(column(i,1),2);
        column(i,3)=xg(column(i,1),1);
        col(j,2)=-1e32;
    end
    
    ll=length(tt);
    sll=zeros(ll,1);
    for i=1:ll
        for j=1:steps
            if tp(j)>tt(i)
                sll(i)=j;
                break;
            end
        end
    end
    
    st0=sll(1);
    for i=1:k
        s0(i,1)=-S((column(i)-1)*4+2,st0);
    end
    
    figure;
    plot(s0,column(:,2),'DisplayName','\sigma_y_0');
    hold on
    pl=zeros(k,1);
    for i=1:ll
        for j=1:k
            pl(j)=Pw(column(j,1),sll(i))-Pw(column(j,1),st0);
        end
        plot(pl,column(:,2),'DisplayName',strcat(num2str(tt(i)),' s'));
    end
    hold off
    legend
    

end

function plot_col(str,PLT)

    load(str,'-mat','GLOBAL','GEOMETRY','SOLVER')

    rr=PLT.range;
    range1=rr(1);
    range2=rr(2);
    
    P=PLT.pw;
    
    if strcmp(PLT.time,'MAX')
        tt=-1;
    else
        tt=PLT.time;
    end
    
    nn=PLT.onnodes;
    if nn==1 && SOLVER.UW==1
        error('only on nodes if the UPw formulation is run')
    end
        

    df=GEOMETRY.df;
    if nn
        x=GEOMETRY.x_0;
        lim=GEOMETRY.nodes;
    else
        x=GEOMETRY.xg_0;
        lim=GEOMETRY.mat_points;
    end

    steps=GLOBAL.ste_p;
    tp=GLOBAL.tp;

    k=0;
    for i=1:lim
        if x(i,1)>range1 && x(i,1)<range2
            k=k+1;
            list_aux(k)=i;
            x_aux(k,:)=x(i,:);
        end
    end

    if k==0
        error('Nodes/mps out of defined range')
    end

    x_col(k,2)=0;
    list_col(k,1)=0;
    for i=1:k
        [~,a]=min(x_aux(:,2));
        x_col(i,:)=x_aux(a,:);
        list_col(i)=a;
        x_aux(a,2)=1e32; 
    end

    if tt~=-1
        ll=length(tt);
        plots=zeros(ll,1);
        for i=1:ll
            for j=1:steps
                if tp(j)>tt(i)
                    plots(i)=j;
                    break;
                end
            end
        end


        figure;
        hold on
        for i=1:length(plots)
            pl=zeros(k,1);
            for j=1:k
                if nn
                    pl(j)=GLOBAL.d((list_col(j)-1)*df+3,plots(i));
                else
                    pl(j)=GLOBAL.pw(list_col(j),plots(i));
                end
            end
            plot(pl,x_col(:,2));
        end

    else
        pl=zeros(k,1);
        for i=1:k
            b=-1e30;
            for j=1:steps
                if nn
                    pp=GLOBAL.d((list_col(i)-1)*df+3,j)/P;
                else
                    pp=GLOBAL.pw(list_col(i),j)/P;
                end
                if pp>b
                    b=pp;
                end
            end
            pl(i)=b;
        end
        figure
        plot(pl,x_col(:,2))
        axis([min(0,min(pl)),max(1.1,max(pl)),...
            min(min(x_col(:,2)),0),max(10,max(x_col(:,2)))])

    end


end

function hh = quadplot(quad,varargin)
%TRIPLOT Plots a 2D triangulation
%   QUADPLOT(QUAD,X,Y) displays the quadrilaterals defined in the
%   M-by-4 matrix QUAD.  A row of QUAD contains indices into X,Y that
%   define a single quadrilateal. The default line color is blue.
%
%   QUADPLOT(...,COLOR) uses the string COLOR as the line color.
%
%   H = QUADPLOT(...) returns a vector of handles to the displayed 
%   quadrilaterals
%
%   QUADPLOT(...,'param','value','param','value'...) allows additional
%   line param/value pairs to be used when creating the plot.
%
%   See also TRISURF, TRIMESH, DELAUNAY, TriRep, DelaunayTri.
%
%   Script code based on copyrighted code from mathworks for TRIPLOT.
%   Allan P. Engsig-Karup, apek@imm.dtu.dk.

error(nargchk(1,inf,nargin,'struct'));

start = 1;

x = varargin{1};
y = varargin{2};
quads = quad;
if (nargin == 3) || (mod(nargin-3,2) == 0)
    c = 'blue';
    start = 3;
else
    c = varargin{3};
    start = 4;
end
  
d = quads(:,[1 2 3 4 1])';
h = plot(x(d), y(d),c,varargin{start:end});
if nargout == 1, hh = h; end
end
