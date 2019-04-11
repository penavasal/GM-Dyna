function PLOT(driver,PLT,varargin)

    folder = PLT.folder;
    freq   = PLT.freq;
    ampl   = PLT.ampl;
    steps  = PLT.steps;
    str_out= PLT.str_out;
    str_in = PLT.str_in;


    if strcmp(driver,'VTK')
        vtk_driver(str_in,str_out,steps,ampl,freq,folder)
    elseif strcmp(driver,'GRID')
        plot_ep(str_in,folder,steps,freq,varargin)
    elseif strcmp(driver,'DEFORMED') || strcmp(driver,'DISTRIBUTION')
        plot_f(driver,str_in,folder,steps,ampl,freq,varargin)
    end
    
end

   
function plot_ep(str,folder,ste_p,rt,varargin)

    
    root=strcat(folder,'/',str,'.mat');
    load(root,'GEOMETRY','SOLVER','GLOBAL');
    
    if strcmp(ste_p,'FULL')
        ste_p=SOLVER.dim;
    end
    
    %[elements,NNE]=size(GEOMETRY.elem);
    x_0=GEOMETRY.x_0;
    
    if strcmp(varargin{1},'GAMMA')
        VAR=GLOBAL.gamma_nds;
    %elseif strcmp(varargin{1},'PW')
    else
        disp('Unrecognized variable')
        stop
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
        pause
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

function plot_f(driver,str,folder,ste_p,hh,rr,varargin)

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
        if strcmp(varargin{1},'PS')
            VAR=GLOBAL.Ps;
        elseif strcmp(varargin{1},'PW')
            VAR=GLOBAL.Pw;
        elseif strcmp(varargin{1},'GAMMA')
            VAR=GLOBAL.gamma;
        else
            disp('Unrecognized variable')
            stop
        end
        if nargin==8
            if strcmp(varargin{2},'FILM')
                film2=1;
            end
        end  
    elseif strcmp(driver,'DEFORMED')
        flag2=1;
        if nargin==7
            if strcmp(varargin{1},'FILM')
                film=1;
            end
        end 
    end

    % Steps
    if strcmp(ste_p,'FULL')
        ste_p=SOLVER.dim;
    end  
    contador=ste_p;
    
    
    % Allocate

    elem=GEOMETRY.elem;
    [Elements,NNE]=size(GEOMETRY.elem);
     x_0=GEOMETRY.x_0;
    [Nodos,~]=size(x_0);
    df=GEOMETRY.df;
    x=zeros(Nodos,1);
    y=zeros(Nodos,1);

    for nd=1:Nodos
        x(nd,1)=x_0(nd,1);
        y(nd,1)=x_0(nd,2);
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
            for nodo=1:Nodos
                if SOLVER.UW==1
                    x2(nodo)=x(nodo)+hh*GLOBAL.d(df*nodo-3,cont);
                    y2(nodo)=y(nodo)+hh*GLOBAL.d(df*nodo-2,cont);
                else
                    x2(nodo)=x(nodo)+hh*GLOBAL.d(df*nodo-1,cont);
                    y2(nodo)=y(nodo)+hh*GLOBAL.d(df*nodo,cont);
                end
            end
            if NNE==3
                triplot(elem,x2,y2)
            elseif NNE==4
                quadplot(elem,x2,y2)
            end
            axis([-0.05,DDD*1.2,-0.05,HHH*1.1])
            drawnow

            if film
                   frame = getframe;
                   writeVideo(movie,frame);  
            end

            if cont==42
                cont
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
                for i=1:Elements
                    for j=1:NNE
                        if nodo==elem(i,j)
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
            for nodo=1:Nodos
                if SOLVER.UW==1
                    x1(nodo)=x(nodo)+hh*GLOBAL.d(df*nodo-3,cont);
                    y1(nodo)=y(nodo)+hh*GLOBAL.d(df*nodo-2,cont);
                else
                    x1(nodo)=x(nodo)+hh*GLOBAL.d(df*nodo-1,cont);
                    y1(nodo)=y(nodo)+hh*GLOBAL.d(df*nodo,cont);
                end
            end
            pr=TriRep(elem,x1(:),y1(:),PWnodo(:,cont));
            trisurf(pr)
            colorbar
            axis([-0.05,DDD*1.2,-0.05,HHH*1.1])
            drawnow

            if film2
                   frame = getframe;
                   writeVideo(movie2,frame);  
            end
            cont
        end

        if film2
            close(movie2);
        end
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