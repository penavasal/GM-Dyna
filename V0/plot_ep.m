    
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