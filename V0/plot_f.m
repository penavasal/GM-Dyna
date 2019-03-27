function plot_f(elem,x_0,ste_p,d,v,Ps,Area,t,Pw,DOF,hh,rr)

film=0; %film or not?
film2=0;
flag2=0; %Displacements of the whole
flag3=1; %Pressure
%hh=1;     %Amplification in pressure plot
%rr=1;   %Frequency in pressure plot

%%            
contador=ste_p;

[Elements,NNE]=size(elem);
[Nodos,sp]=size(x_0);
df=DOF*sp;
x=zeros(Nodos,1);
y=zeros(Nodos,1);

for nd=1:Nodos
    x(nd,1)=x_0(nd,1);
    y(nd,1)=x_0(nd,2);
end

DDD=max(x(:));
HHH=max(y(:));


if film
    movie = VideoWriter('video.avi');
    open(movie);
end
if film2
    movie2 = VideoWriter('press.avi');
    open(movie2);
end


%%
if flag2==1
    x2=zeros(Nodos,1);
    y2=zeros(Nodos,1);
    figure
    for cont=1:rr:contador
        for nodo=1:Nodos
            if DOF==2
                x2(nodo)=x(nodo)+hh*d(df*nodo-3,cont);
                y2(nodo)=y(nodo)+hh*d(df*nodo-2,cont);
            else
                x2(nodo)=x(nodo)+hh*d(df*nodo-1,cont);
                y2(nodo)=y(nodo)+hh*d(df*nodo,cont);
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

%% Pressure evolution 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag3
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
                        suma1=suma1+Pw(i,cont)*Area(i);
                        suma2=suma2+Area(i);
                    end
                end
            end;
            PWnodo(nodo,cont)=suma1/suma2;
        end;
    end;
    
    
    figure
    for cont=1:rr:ste_p
        for i=1:Nodos
            x1(i)=x(i)+hh*d(i*df-3,cont);
            y1(i)=y(i)+hh*d(i*df-2,cont);
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
    end;
    
    if film2
        close(movie2);
    end
end


end
