function read_geometry(UW,ELEMENT,DIM,PLOT_ini,AMP,filename,pathgeo)

    global GEOMETRY
    
    [x,elem,NNE,El_type]=read_mesh(filename,pathgeo,DIM);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial checks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(ELEMENT{1},'L1')
        if DIM>1
            fprintf('Error, DIM mismatches with spatial dimension of the mesh!!\n')
            stop
        end
    else
        if strcmp(ELEMENT{1},'L6')
            f1=strcat(filename,'_esq');
            if ismac || isunix  % Code to run on Mac or Unix plaform 
                f_i=strcat(pathgeo,'/',f1,'.msh');
            elseif ispc         % Code to run on Windows platform
                f_i=strcat(pathgeo,'\',f1,'.msh');
            else
                disp('Platform not supported')
                stop
            end 
            if isfile(f_i)
                [x_e,elemesq,~,~]=Geo_DynCLM(filename,pathgeo,DIM);
            else
                PLOT_ini=0;
            end
        end
            
        if DIM==3
            fprintf('Error, DIM mismatches with spatial dimension of the mesh!!\n')
            stop
        end 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build the mesh
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_a=x*AMP;     
    if strcmp(ELEMENT{1},'T3') ...
            || strcmp(ELEMENT{1},'T3-inverse')...
            || strcmp(ELEMENT{1},'T3-diamond')
        
        if NNE==4
            if strcmp(ELEMENT{1},'T3')
                [GEOMETRY.elem,GEOMETRY.patch_con,GEOMETRY.patch_el]=...
                    split2(x_a,elem);
                [xg,GEOMETRY.Area]=g_center(x_a,GEOMETRY.elem,DIM);
                elem_0=GEOMETRY.elem;
            elseif strcmp(ELEMENT{1},'T3-inverse')
                [elem]=reverse(elem); %reverse
                [GEOMETRY.elem,GEOMETRY.patch_con,GEOMETRY.patch_el]=...
                    split2(x_a,elem);
                [xg,GEOMETRY.Area]=g_center(x_a,GEOMETRY.elem,DIM);
                elem_0=GEOMETRY.elem;
            elseif strcmp(ELEMENT{1},'T3-diamond')
                %Patch 4P1P0
                [xp]=q_g_center(x_a,elem);
                [x_a,GEOMETRY.elem,GEOMETRY.patch_con,GEOMETRY.patch_el]=...
                    split(x_a,elem,xp);
                [xg,GEOMETRY.Area]=g_center(x_a,GEOMETRY.elem,DIM);
                elem_0=GEOMETRY.elem;
            end
            
        elseif NNE==3
            GEOMETRY.elem=elem;
            [GEOMETRY.ELEMENT(:),GEOMETRY.ELEMENT(:).Area]=...
                g_center(x_a,GEOMETRY.elem,DIM);
            [elements,~]=size(GEOMETRY.elem);
            GEOMETRY.x_0=x_a;
            elem_0=GEOMETRY.elem;
            GEOMETRY.patch_el=(1:elements)';  
            GEOMETRY.patch_con=(1:elements)'; 
        end
        
    elseif strcmp(ELEMENT{1},'Q4')
        
        GEOMETRY.elem=elem;
        [elements,~]=size(elem);
        [xg,GEOMETRY.Area]=g_center(x_a,elem,DIM);
        elem_0=GEOMETRY.elem;
        %No patches
        GEOMETRY.patch_el=(1:elements)';  
        GEOMETRY.patch_con=(1:elements)'; 
        
    elseif strcmp(ELEMENT{1},'Q4-4')
        elem_0=elem;
        [xg,GEOMETRY.elem,GEOMETRY.Area,GEOMETRY.patch_el,...
            GEOMETRY.patch_con]=quad4xg(x_a,elem);
        
    elseif strcmp(ELEMENT{1},'T6')
        [xg,GEOMETRY.Area]=g_center(x_a,elemesq,DIM);
        [elements,~]=size(GEOMETRY.elem);
        GEOMETRY.patch_el=(1:elements)';  
        GEOMETRY.patch_con=(1:elements)';
        
    end

    [aa,bb]=size(GEOMETRY.patch_con);
    GEOMETRY.Area_p=zeros(aa,1);
    for i=1:aa
        for j=1:bb
            GEOMETRY.Area_p(i)=GEOMETRY.Area_p(i) +...
                GEOMETRY.Area(GEOMETRY.patch_con(i,j));
        end
    end 
    
    GEOMETRY.x_0=x_a;
    GEOMETRY.xg=xg;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Geometry parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [GEOMETRY.nodes,GEOMETRY.sp]=size(GEOMETRY.x_0);
    [GEOMETRY.elements,NNE_f]=size(GEOMETRY.elem);
    
    % Degrees of freedom
    if UW==0
        GEOMETRY.df=GEOMETRY.sp;
    elseif UW==1
        GEOMETRY.df=2*GEOMETRY.sp;
    elseif UW==2
        GEOMETRY.df=GEOMETRY.sp+1;
    elseif UW==3
        GEOMETRY.df=2*GEOMETRY.sp+1;
    end   
    
    % Mesh size
    h(GEOMETRY.elements,1)=0;
    GEOMETRY.h_nds(GEOMETRY.nodes,1)=0;
    for e=1:GEOMETRY.elements
        if NNE_f==3
            h(e)=sqrt(2*GEOMETRY.Area(e));
        elseif NNE_f==4
            h(e)=sqrt(GEOMETRY.Area(e));
        end
        for i=1:NNE_f
            GEOMETRY.h_nds(GEOMETRY.elem(e,i))=...
                max(GEOMETRY.h_nds(GEOMETRY.elem(e,i)),h(e));
        end
    end
    GEOMETRY.h_ini=h;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if PLOT_ini
        plot_nb(0,0,x_a,xg,elem_0,0,0)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add shape function parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(ELEMENT{2},'LME')
        read_LME;
    end

end

function [x_a,nw_elem,patch_con,patch_el]=split(x_a,elem,xp)

    [nodes,sp]=size(x_a);
    [elements,NNE]=size(elem);
    
    for i=1:elements
        for j=1:sp
            x_a(nodes+i,j)=xp(i,j);
        end
    end
        
    nw_elem=zeros(elements*NNE,NNE-1);
    patch_el=zeros(elements,1);
    patch_con=zeros(elements,NNE);
    for i=1:elements
        for j=1:NNE-1
            nw_elem((i-1)*NNE+j,1)=nodes+i;
            nw_elem((i-1)*NNE+j,2)=elem(i,j);
            nw_elem((i-1)*NNE+j,3)=elem(i,j+1);
        end
        j=j+1;
        nw_elem((i-1)*NNE+j,1)=nodes+i;
        nw_elem((i-1)*NNE+j,2)=elem(i,j);
        nw_elem((i-1)*NNE+j,3)=elem(i,1);
        
        for j=1:NNE
            patch_con(i,j)=(i-1)*NNE+j;
        end
    end
    
    for j=1:elements
        for k=1:NNE
            patch_el(patch_con(j,k))=j;
        end
    end
end

function [nw_elem,patch_con,patch_el]=split2(x_a,elem)

    [elements,NNE]=size(elem);
    [~,sp]=size(x_a);
        
    nw_elem=zeros(elements*sp,NNE-1);
    patch_el=zeros(elements,1);
    patch_con=zeros(elements,sp);

    for i=1:elements
        X=zeros(NNE,sp);
        for j=1:NNE
            for k=1:sp
                X(j,k)=x_a(elem(i,j),k);
            end
        end
        [t]=circunf3(X);
        
        if t==1
            nw_elem((i-1)*sp+1,1)=elem(i,1);
            nw_elem((i-1)*sp+1,2)=elem(i,2);
            nw_elem((i-1)*sp+1,3)=elem(i,3);

            nw_elem((i-1)*sp+2,1)=elem(i,3);
            nw_elem((i-1)*sp+2,2)=elem(i,4);
            nw_elem((i-1)*sp+2,3)=elem(i,1);
        else
            nw_elem((i-1)*sp+1,1)=elem(i,2);
            nw_elem((i-1)*sp+1,2)=elem(i,3);
            nw_elem((i-1)*sp+1,3)=elem(i,4);

            nw_elem((i-1)*sp+2,1)=elem(i,4);
            nw_elem((i-1)*sp+2,2)=elem(i,1);
            nw_elem((i-1)*sp+2,3)=elem(i,2);
        end
        
        for j=1:sp
            patch_con(i,j)=(i-1)*sp+j;
        end
    end
    
    for j=1:elements
        for k=1:sp
            patch_el(patch_con(j,k))=j;
        end
    end
end

function [t]=circunf3(x_a)
 
    A=ones(3);
    C=zeros(3,1);
    t=0;
    
    for i=1:3
        C(i)=-(x_a(i,1)^2+x_a(i,2)^2);
        for j=2:1:3
            A(i,j)=x_a(i,j-1);
        end
    end
    
    B=A\C;
    
    c=[-B(2)/2; -B(3)/2];
    r=sqrt(c(1)^2+c(2)^2-B(1));
    
    
    d=sqrt((x_a(4,1)-c(1))^2+(x_a(4,2)-c(2))^2);
    if d>=r-r*0.05
        t=1;
    end
    
end

function [elem2]=reverse(elem)

    [elements,NNE]=size(elem);
    elem2=zeros(elements,NNE);
    for i=1:elements
        elem2(i,1)=elem(i,4);
        elem2(i,2)=elem(i,1);
        elem2(i,3)=elem(i,2);
        elem2(i,4)=elem(i,3);
    end

end

function [xg,elem,Area,patch_el,patch_con]=quad4xg(x_a,elem_0)

    [~,sp]=size(x_a);
    [elements,NNE]=size(elem_0);
    
    [x0]=q_g_center(x_a,elem_0);
    
    xg=zeros(elements*NNE,sp);
    Area=zeros(elements*NNE,1);
    patch_el=zeros(elements*NNE,1);
    patch_con=zeros(elements,NNE);
    k=0;
    for e=1:elements
        xn=zeros(NNE+1,1);
        yn=zeros(NNE+1,1);

        for i=1:NNE
            nd=elem_0(e,i);
            xn(i)=x_a(nd,1);
            yn(i)=x_a(nd,2);
        end
        nd=elem_0(e,1);
        xn(i+1)=x_a(nd,1);
        yn(i+1)=x_a(nd,2);

        a1 = 0;
        a2 = 0;

        for i=1:NNE
            a1 = a1 + xn(i)*yn(i+1);
            a2 = a2 + yn(i)*xn(i+1);
        end
        area_=abs(a1-a2)/2;
        
        for i=1:NNE
            xb=x_a(elem_0(e,i),1);
            yb=x_a(elem_0(e,i),2);
            
            MR=(yb-x0(e,2))/(xb-x0(e,1));
            MS=-(xb-x0(e,1))/(yb-x0(e,2));
            
            %Mat=inv([-MS 1; -MR 1]);
            Mat=[-MS 1; -MR 1];
            
            d=sqrt( (xb-x0(e,1))^2 + (yb-x0(e,2))^2 )/sqrt(3);
            
            if (yb-x0(e,2))<0
                f=[ x0(e,2) - MS*x0(e,1) - sqrt(MS^2+1)*d;
                    x0(e,2) - MR*x0(e,1)];
            else  
                f=[ x0(e,2) - MS*x0(e,1) + sqrt(MS^2+1)*d;
                    x0(e,2) - MR*x0(e,1)];
            end
            xp=Mat\f;
            
            %f2=[ x0(e,2) - MS*x0(e,1) + sqrt(MS^2+1)*d;
            %    x0(e,2) - MR*x0(e,1)];
            %xp2=Mat*f2;
            
            k=k+1;
            
            for j=1:sp
                xg(k,j)=xp(j);
            end
            elem(k,:)=elem_0(e,:);
            Area(k)=area_/NNE;
            
        end 
        for j=1:4
            patch_con(e,j)=(e-1)*NNE+j;
        end
    end
    
    for j=1:elements
        for k=1:NNE
            patch_el(patch_con(j,k))=j;
        end
    end
end

function [xg]=q_g_center(x_a,elem)

    [elements,NNE]=size(elem);
    [~,sp]=size(x_a);

    xg=zeros(elements,sp);
    
    for e=1:elements
        sum=zeros(sp);
        for i=1:NNE
            nd=elem(e,i);
            for j=1:sp                    
                sum(j)=sum(j)+x_a(nd,j);
            end
        end
        for j=1:sp
            xg(e,j)=sum(j)/NNE;
        end
    end
end

function [xg,Area]=g_center(x_a,elem,D1)

    [elements,NNE]=size(elem);
    [~,sp]=size(x_a);

    xg=zeros(elements,sp);
    Area=zeros(elements,1);
    
    for e=1:elements
        sum=zeros(sp);
        for i=1:NNE
            nd=elem(e,i);
            for j=1:sp                    
                sum(j)=sum(j)+x_a(nd,j);
            end
        end
        for j=1:sp
            xg(e,j)=sum(j)/NNE;
        end
    end
    
   if D1~=1
        for e=1:elements
            xn=zeros(NNE+1,1);
            yn=zeros(NNE+1,1);

            for i=1:NNE
                nd=elem(e,i);
                xn(i)=x_a(nd,1);
                yn(i)=x_a(nd,2);
            end
            nd=elem(e,1);
            xn(i+1)=x_a(nd,1);
            yn(i+1)=x_a(nd,2);

            a1 = 0;
            a2 = 0;

            for i=1:NNE
                a1 = a1 + xn(i)*yn(i+1);
                a2 = a2 + yn(i)*xn(i+1);
            end
            Area(e)=abs(a1-a2)/2;
        end
    else
        for e=1:elements
        	Area(e)=abs(x_a(elem(e,1),1)-x_a(elem(e,2),1));
        end
   end
        
end

function [x,elem,NNE,El_type]=read_mesh(str1,str2,DIM)


    close all
    
    % File: Geo_DynCLM
    %   Read geometry in *.msh and save in the geometry variables

    % Date:
    %   Version 2.0   02.06.2019

    %Input file name
    if ismac || isunix  % Code to run on Mac or Unix plaform 
        filename_i=strcat(str2,'/',str1,'.msh');
    elseif ispc         % Code to run on Windows platform
        filename_i=strcat(str2,'\',str1,'.msh');
    else
        disp('Platform not supported')
        stop
    end 

    f_i = fopen(filename_i, 'rt');
    formato = '%s %s %s %s %s %s %s'; % formato de cada línea 

    tsargs = {...
        'HeaderLines',0,...
        'HeaderColumns',0,...
        'ReturnOnError',false,...
        'EmptyValue',0,...
        'CollectOutput',true,...
        'EndOfLine','\r\n'};

    % Collect info
    res  = textscan(f_i,formato,-1,tsargs{:});
    data=res{1};
    
    fclose(f_i);
    
    % CHECK
    if DIM>str2double(data{1,3})
    	fprintf('Error, DIM mismatches with spatial dimension of the mesh!!\n')
        stop
    end
    
    El_type=data{1,5};
    NNE=str2double(data{1,7});
        

    fin=0;
    t=3;
    in=1;
    while fin==0
        if strcmp(data{t,1},'End')
            break;
        end
        for i=2:DIM+1
            x(in,i-1)=str2double(data{t,i});
        end
        t=t+1;
        in=in+1;
    end
    
    t=t+2;
    in=1;
    while fin==0
        if strcmp(data{t,1},'End')
            break;
        end
        for i=2:NNE+1
            elem(in,i-1)=str2double(data{t,i});
        end
        t=t+1;
        in=in+1;
    end

end