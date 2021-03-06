

function plot_nb(NODO,nds,x_0,elem,ds,h)

global GEOMETRY

xg=GEOMETRY.xg_0;

[~,NNE]=size(elem);
[nodes,sp]=size(x_0);
x_a=x_0;
if h~=0
    for i=1:nodes
        x_a(i,1)=x_a(i,1)+h*ds(i*sp-1,1);
        x_a(i,2)=x_a(i,2)+h*ds(i*sp,1);
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
