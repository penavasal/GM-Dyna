function [MAT_POINT,NODE_LIST]=...
    read_geometry(ELEMENT,GRID,DIM,PLOT_ini,AMP,filename,filegrid,pathgeo)

    global GEOMETRY SOLVER
    
    UW=SOLVER.UW;
    AXI=SOLVER.AXI;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRID TYPE
    % Explore if there is a different grid or we employ the original file
    % Initial checks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    grid=0;
    
    if DIM==1
        if not(strcmp(GRID,'L1'))
            fprintf('Error, DIM mismatches with spatial dimension of the mesh!!\n')
            stop
        end
    elseif DIM>2
            fprintf('Error, DIM not implemented!!\n')
            stop
    else
        if not(strcmp(filegrid,filename))
            if isempty(filegrid)
                filegrid=filename;
                if isempty(GRID) 
                    if SOLVER.TYPE{1}==1
                        fprintf('Error, Undefined grid type!!\n')
                        stop
                    else
                        if strcmp(ELEMENT,'T3') || strcmp(ELEMENT,'T3-inverse')...
                                || strcmp(ELEMENT,'T3-diamond') 
                            GRID='T3';
                        elseif strcmp(ELEMENT,'T3-3')
                            GRID='T3';
                        elseif strcmp(ELEMENT,'Q4') || strcmp(ELEMENT,'Q4-4')
                            GRID='Q4';
                        elseif strcmp(ELEMENT,'T6') || strcmp(ELEMENT,'T6-3')
                            GRID='T6';
                        elseif strcmp(ELEMENT,'T6P3') || strcmp(ELEMENT,'T6P3-3')
                            if UW==0
                                disp('water-type element, cannot be U formulation');
                                stop;
                            else
                                GRID='T6';
                            end
                        elseif strcmp(ELEMENT,'Q8') || strcmp(ELEMENT,'Q8-4')
                            GRID='Q8';
                        elseif strcmp(ELEMENT,'Q8P4') || strcmp(ELEMENT,'Q8P4-4')
                            if UW==0
                                disp('water-type element, cannot be U formulation');
                                stop;
                            else
                                GRID='Q8';
                            end
                        elseif strcmp(ELEMENT,'L1')
                            GRID='L1';
                        else
                        	fprintf('Error, Undefined grid type!!\n')
                            stop
                        end
                    end
                elseif strcmp(GRID,'Q4')
                    if strcmp(ELEMENT,'T3') || strcmp(ELEMENT,'L1') ||...
                    strcmp(ELEMENT,'T3-inverse')|| strcmp(ELEMENT,'T3-diamond')
                        fprintf('Error, wrong definition of grid type!!\n')
                        stop
                    end
                elseif strcmp(GRID,'T6') || strcmp(GRID,'T3')
                    if strcmp(ELEMENT,'Q4') || strcmp(ELEMENT,'Q4-4') ||...
                            strcmp(ELEMENT,'L1')
                        fprintf('Error, wrong definition of grid type!!\n')
                        stop
                    end
                elseif strcmp(GRID,'L1')
                    if not(strcmp(ELEMENT,'L1'))
                        fprintf('Error, wrong definition of grid type!!\n')
                        stop
                    end  
                end
            elseif SOLVER.TYPE{1}~=1
                filegrid=filename;
            else
                grid=1;
                [x_mp,elem_mp,NNE_mp,~]=read_mesh(filename,pathgeo,DIM);
            end
        else
            if strcmp(GRID,'T6')
                f1=strcat(filegrid,'_esq');
                if ismac || isunix  % Code to run on Mac or Unix plaform 
                    f_i=strcat(pathgeo,'/',f1,'.msh');
                elseif ispc         % Code to run on Windows platform
                    f_i=strcat(pathgeo,'\',f1,'.msh');
                else
                    disp('Platform not supported')
                    stop
                end 
                if isfile(f_i)
                    [x_e,elem_0,~,~]=read_mesh(f_i,pathgeo,DIM);
                else
                    disp('Cannot plot initial mesh, check it please')
                    PLOT_ini=0;
                end
            end
        end
            
    end


    GEOMETRY=struct('elem',0,'elem_c',0,'patch_con',0,'patch_el',0,...
        'Area',0,'x_0',0,'ELEMENT',ELEMENT,...
        'elements',0,'nodes',0,'Area_p',0,'sp',0,'df',0,'xg_0',0,...
        'h_ini',0,'h_nds',0,'mat_points',0,'node_connect',0,...
        'material',0);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build the grid mesh
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
    
    
    %[x,elem,NNE,~]=read_mesh(filegrid,pathgeo,DIM);
    [x,elem,NNE,materials,NODE_LIST]=read_dat(filegrid,pathgeo,DIM);
    x_a=x*AMP;    % SCALE
    
    if strcmp(GRID,'T3')
        
        if NNE==4 && grid==0
            if strcmp(ELEMENT,'T3-3')
                [xg,GEOMETRY.Area,elem,patch_con,patch_el,materials]=...
                    split3xg(x_a,elem,materials);
            elseif strcmp(ELEMENT,'T3')
                [elem,patch_con,patch_el,materials]=...
                    split2(x_a,elem,materials);
            elseif strcmp(ELEMENT,'T3-inverse')
                [elem]=reverse(elem); %reverse
                [elem,patch_con,patch_el,materials]=...
                    split2(x_a,elem,materials);
            elseif strcmp(ELEMENT,'T3-diamond')
                %Patch 4P1P0
                [xp]=q_g_center(x_a,elem);
                [x_a,elem,patch_con,patch_el,materials]=...
                    split(x_a,elem,xp,materials);
            end
            
            GEOMETRY.elem=elem;
            GEOMETRY.elem_c=elem;
            GEOMETRY.patch_con=patch_con;
            GEOMETRY.patch_el=patch_el;
            if strcmp(ELEMENT,'T3-3')
                
            else
                [xg,GEOMETRY.Area]=g_center(x_a,GEOMETRY.elem,DIM);
            end
            
        elseif NNE==3
            GEOMETRY.elem=elem;
            GEOMETRY.elem_c=elem;
            if grid==0
                if strcmp(ELEMENT,'T3-3')
                    [xg,GEOMETRY.Area,GEOMETRY.patch_el,...
                        GEOMETRY.patch_con,materials]=tri3xg(x_a,elem,materials);
                else
                    [xg,GEOMETRY.Area]=g_center(x_a,GEOMETRY.elem,DIM);
                     [elements,~]=size(GEOMETRY.elem);
                    GEOMETRY.patch_el=(1:elements)';  
                    GEOMETRY.patch_con=(1:elements)'; 
                end
                GEOMETRY.x_0=x_a;
            end
        end
            
    elseif strcmp(GRID,'T6')
        if grid==1
            GEOMETRY.elem=elem;
            GEOMETRY.elem_c=elem;
        else 
            if strcmp(ELEMENT,'T6') || strcmp(ELEMENT,'T6P3')
                if NNE~=6 
                    disp('NNEnot equal to 6, incorrect');
                    stop;
                end
                [elem]=corner_nds6(elem);
                elem_c=elem(:,1:3);
                GEOMETRY.elem_c=elem_c;
                GEOMETRY.elem=elem;
                [xg,GEOMETRY.Area]=g_center(x_a,elem_c,DIM);
                [elements,~]=size(GEOMETRY.elem);
                GEOMETRY.patch_el=(1:elements)';  
                GEOMETRY.patch_con=(1:elements)';
            elseif strcmp(ELEMENT,'T6-3') || strcmp(ELEMENT,'T6P3-3')
                if NNE~=6 
                    disp('NNEnot equal to 6, incorrect');
                    stop;
                end
                [elem]=corner_nds6(elem);
                elem_c=elem(:,1:3);
                GEOMETRY.elem_c=elem_c;
                GEOMETRY.elem=elem;
                [xg,GEOMETRY.Area,GEOMETRY.patch_el,...
                    GEOMETRY.patch_con,materials]=tri3xg(x_a,elem_c,materials);
            end
        end
        
    elseif strcmp(GRID,'Q4') || strcmp(GRID,'Q8')
        if grid==1
            GEOMETRY.elem=elem;
            GEOMETRY.elem_c=elem;
        else
            if strcmp(ELEMENT,'Q4') || strcmp(ELEMENT,'Q8') || strcmp(ELEMENT,'Q8P4')
                [elements,NNE]=size(elem);
                if NNE==8 && strcmp(ELEMENT,'Q4')
                    disp('NNE=8 and element Q4, incorrect');
                    stop;
                end
                if strcmp(ELEMENT,'Q8') || strcmp(ELEMENT,'Q8P4')
                    [elem]=corner_nds(elem);
                    elem_c=elem(:,1:4);
                else
                    elem_c=elem;
                end
                GEOMETRY.elem_c=elem_c;
                [xg,GEOMETRY.Area]=g_center(x_a,elem_c,DIM);
                GEOMETRY.patch_el=(1:elements)';  
                GEOMETRY.patch_con=(1:elements)'; 
                GEOMETRY.elem=elem;
            elseif strcmp(ELEMENT,'Q4-4') || strcmp(ELEMENT,'Q8-4') || strcmp(ELEMENT,'Q8P4-4')
                if NNE==8 && strcmp(ELEMENT,'Q4-4')
                    disp('NNE=8 and element Q4-4, incorrect');
                    stop;
                end
                if strcmp(ELEMENT,'Q8-4') || strcmp(ELEMENT,'Q8P4-4')
                    [elem]=corner_nds(elem);
                    elem_c=elem(:,1:4);
                else
                    elem_c=elem;  
                end
                GEOMETRY.elem_c=elem_c;
                GEOMETRY.elem=elem;
                [xg,GEOMETRY.Area,GEOMETRY.patch_el,...
                    GEOMETRY.patch_con,materials]=quad4xg(x_a,elem_c,materials);
            end
        end
    end
    
    if not(strcmp(ELEMENT,'T6'))
        x_e=x_a;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get the material point position
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if grid==1
        x_mp=x_mp*AMP;     
        if strcmp(ELEMENT,'T3') ||  strcmp(ELEMENT,'T3-inverse') ||...
                strcmp(ELEMENT,'T3-diamond')
            if NNE_mp==4
                if strcmp(ELEMENT,'T3')
                    [~,GEOMETRY.patch_con,GEOMETRY.patch_el]=...
                        split2(x_mp,elem_mp);
                elseif strcmp(ELEMENT,'T3-inverse')
                    [elem_mp]=reverse(elem_mp); %reverse
                    [~,GEOMETRY.patch_con,GEOMETRY.patch_el]=...
                        split2(x_mp,elem_mp);
                elseif strcmp(ELEMENT,'T3-diamond')
                    %Patch 4P1P0
                    [xp]=q_g_center(x_mp,elem_mp);
                    [x_mp,~,GEOMETRY.patch_con,GEOMETRY.patch_el]=...
                        split(x_mp,elem_mp,xp);
                end
                
                [xg,GEOMETRY.Area]=g_center(x_mp,elem_mp,DIM);

            elseif NNE_mp==3
                [xg,GEOMETRY.Area]=g_center(x_mp,elem_mp,DIM);
                [els,~]=size(elem_mp);
                GEOMETRY.patch_el=(1:els)';  
                GEOMETRY.patch_con=(1:els)'; 
            end

        elseif strcmp(ELEMENT,'Q4')

            [els,~]=size(elem_mp);
            [xg,GEOMETRY.Area]=g_center(x_mp,elem_mp,DIM);
            %No patches
            GEOMETRY.patch_el=(1:els)';  
            GEOMETRY.patch_con=(1:els)'; 

        elseif strcmp(ELEMENT,'Q4-4')
            [xg,GEOMETRY.Area,GEOMETRY.patch_el,GEOMETRY.patch_con,...
                materials]=quad4xg(x_mp,elem_mp,materials);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Geometry parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    GEOMETRY.material=materials;
    [aa,bb]=size(GEOMETRY.patch_con);
    GEOMETRY.Area_p=zeros(aa,1);
    for i=1:aa
        for j=1:bb
            GEOMETRY.Area_p(i)=GEOMETRY.Area_p(i) +...
                GEOMETRY.Area(GEOMETRY.patch_con(i,j));
        end
    end 
    
    
    GEOMETRY.x_0  = x_a;
    GEOMETRY.xg_0 = xg;
    GEOMETRY.H    = max(xg(:,2));
    
    [GEOMETRY.nodes,GEOMETRY.sp]=size(GEOMETRY.x_0);
    [GEOMETRY.elements,NNE_f]=size(GEOMETRY.elem);
    
    [GEOMETRY.mat_points,~]=size(xg);
    
    % Degrees of freedom
    if UW==0
        GEOMETRY.df=GEOMETRY.sp;
    elseif UW==1 || UW==4
        GEOMETRY.df=2*GEOMETRY.sp;
    elseif UW==2
        GEOMETRY.df=GEOMETRY.sp+1;
    elseif UW==3
        GEOMETRY.df=2*GEOMETRY.sp+1;
    end
    
    if AXI==0
        if GEOMETRY.sp==2
            GEOMETRY.b_dim=3;
            GEOMETRY.s_dim=4;
            GEOMETRY.f_dim=5;
        elseif GEOMETRY.sp==1
            GEOMETRY.b_dim=1;
            GEOMETRY.s_dim=1;
            GEOMETRY.f_dim=1;
        else
            GEOMETRY.b_dim=6;
            GEOMETRY.s_dim=6;
            GEOMETRY.f_dim=9;
        end  
    elseif AXI==1
        GEOMETRY.b_dim=4;
        GEOMETRY.s_dim=4;
        GEOMETRY.f_dim=5;
    else 
        disp('DIMENSION OF THE AXI PARAMETER??')
        stop
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creation of MAT_POINT object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    % Dimension of the NNE depending on the support
    
    
    [phases,~]=size(SOLVER.PHASES);
    NNEp=zeros(phases,1);
    for i=1:phases
        NNEp(i)=NNE_f;
        if i==2
            if strcmp(ELEMENT,'Q8P4-4') || strcmp(ELEMENT,'Q8P4')
                NNEp(i)=4;
            elseif strcmp(ELEMENT,'T6P3-3') || strcmp(ELEMENT,'T6P3')
                NNEp(i)=3;
            end
        end
        
        MAT_POINT{i}(1:GEOMETRY.mat_points)=...
                  struct('xg', zeros(GEOMETRY.sp,1),...
                         'element', zeros(1),...
                         'near', zeros(NNEp(i),1),...
                         'N', zeros(NNEp(i),1),...
                         'B', zeros(GEOMETRY.b_dim,GEOMETRY.sp*NNEp(i)),...
                         'EP', zeros(3,2),...
                         'J', ones(1),...
                         'w', zeros(1),...
                         'xi', zeros(GEOMETRY.sp,1),...
                         'near_p',zeros(1),...
                         'REM',zeros(3,1)...
                     );
    end
             
    [MAT_POINT{1}]=LIB.list2S(MAT_POINT{1},'xg',xg);
    
    %-----------------------------------------
    % NODE CONNECTIVITY
    %-----------------------------------------
    node_connect=cell(GEOMETRY.nodes,1);
    for e=1:GEOMETRY.elements
        for i=1:NNE_f
            node_connect(GEOMETRY.elem(e,i))=...
                {cat(2,node_connect{GEOMETRY.elem(e,i)},e)};
        end
    end
    GEOMETRY.node_connect=node_connect;
    
    %--------------------------------------------------
    % near ELEMENTs
    %-------------------------------------------------- 
    element_near=cell(GEOMETRY.elements,1);
    for e=1:GEOMETRY.elements
        for i=1:NNE_f
            nc=node_connect{GEOMETRY.elem(e,i)};
            for j=1:length(nc)
                if nc(j)~=e && not(ismember(nc(j),element_near{e}))
                	element_near(e)=...
                        {cat(2,element_near{e},nc(j))};
                end
            end
        end
        GEOMETRY.element_near{e,1}=sort(element_near{e});
    end
    %--------------------------------------------------
    % MATERIAL POINT in ELEMENT, and initial neighbors
    %--------------------------------------------------       
    if grid==1
        for mp=1:GEOMETRY.mat_points
            if I==0
                e=0;
                while I==0
                    e=e+1;
                    if e>GEOMETRY.elements
                        disp('Error, out of the mesh');
                        stop
                    end
                    [I]=LIB.IoO(getfield(S(mp),'xg'),...
                            GEOMETRY.x_a,GEOMETRY.elem(e,:));  
                end 
            end
            MAT_POINT{1}(mp).element=e;
            MAT_POINT{1}(mp).near=GEOMETRY.elem(e,:);
        end
    else
        if strcmp(ELEMENT,'Q4-4') || strcmp(ELEMENT,'Q8-4') || strcmp(ELEMENT,'Q8P4-4')
            mp=0;
            for e=1:GEOMETRY.elements
                for j=1:4
                    mp=mp+1;
                    MAT_POINT{1}(mp).element=e;
                    MAT_POINT{1}(mp).near=GEOMETRY.elem(e,:);
                end
            end
        elseif strcmp(ELEMENT,'T3-3') || strcmp(ELEMENT,'T6-3') || strcmp(ELEMENT,'T6P3-3')
            mp=0;
            for e=1:GEOMETRY.elements
                for j=1:3
                    mp=mp+1;
                    MAT_POINT{1}(mp).element=e;
                    MAT_POINT{1}(mp).near=GEOMETRY.elem(e,:);
                end
            end
        else
            for mp=1:GEOMETRY.mat_points
                MAT_POINT{1}(mp).element=mp;
                MAT_POINT{1}(mp).near=GEOMETRY.elem(mp,:);
            end
        end
    end
    
    for i=2:phases
        for mp=1:GEOMETRY.mat_points
            MAT_POINT{i}(mp).element=MAT_POINT{1}(mp).element;
            MAT_POINT{i}(mp).xg=MAT_POINT{1}(mp).xg;
            if strcmp(ELEMENT,'Q8P4-4') || strcmp(ELEMENT,'Q8P4')...
                    || strcmp(ELEMENT,'T6P3') || strcmp(ELEMENT,'T6P3-3')
                e=MAT_POINT{1}(mp).element;
                MAT_POINT{i}(mp).near=GEOMETRY.elem_c(e,:);
            else
                MAT_POINT{i}(mp).near=MAT_POINT{1}(mp).near;
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mesh size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    GEOMETRY.h_ini=zeros(GEOMETRY.mat_points,1);
    GEOMETRY.h_nds(GEOMETRY.nodes,1)=0;
    aux=zeros(GEOMETRY.nodes,1);
    for e=1:GEOMETRY.mat_points
        if strcmp(ELEMENT,'T3') || strcmp(ELEMENT,'T3-inverse') || strcmp(ELEMENT,'T3-3')...
                || strcmp(ELEMENT,'T3-diamond') || strcmp(ELEMENT,'T6')...
                ||strcmp(ELEMENT,'T6-3') ||strcmp(ELEMENT,'T6P3-3') ||strcmp(ELEMENT,'T6P3')
            GEOMETRY.h_ini(e)=sqrt(2*GEOMETRY.Area(e));
        elseif strcmp(ELEMENT,'Q4-4')||strcmp(ELEMENT,'Q4') ||strcmp(ELEMENT,'Q8')...
                 || strcmp(ELEMENT,'Q8P4')|| strcmp(ELEMENT,'Q8P4-4')...
                 || strcmp(ELEMENT,'Q8-4')
            GEOMETRY.h_ini(e)=sqrt(GEOMETRY.Area(e));
        end
        near=MAT_POINT{1}(e).near;
        for i=1:length(near)
            aux(near(i))=aux(near(i))+1;
            GEOMETRY.h_nds(near(i))=...
                GEOMETRY.h_nds(near(i))+GEOMETRY.h_ini(e);
        end
    end
    
    GEOMETRY.h_nds=GEOMETRY.h_nds./aux;
    GEOMETRY.Area=GEOMETRY.Area*SOLVER.thickness;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if PLOT_ini
        PLOT('NEIGHBORS',pwd,'AMPL',0);
    end


end

function [x_a,nw_elem,patch_con,patch_el,nw_mat]=split(x_a,elem,xp,mat)

    [nodes,sp]=size(x_a);
    [elements,NNE]=size(elem);
    
    for i=1:elements
        for j=1:sp
            x_a(nodes+i,j)=xp(i,j);
        end
    end
        
    nw_elem=zeros(elements*NNE,NNE-1);
    nw_mat=zeros(elements*NNE,1);
    patch_el=zeros(elements,1);
    patch_con=zeros(elements,NNE);
    for i=1:elements
        for j=1:NNE-1
            nw_elem((i-1)*NNE+j,1)=nodes+i;
            nw_elem((i-1)*NNE+j,2)=elem(i,j);
            nw_elem((i-1)*NNE+j,3)=elem(i,j+1);
            
            nw_mat((i-1)*NNE+j,1)=mat(i);
        end
        j=j+1;
        nw_elem((i-1)*NNE+j,1)=nodes+i;
        nw_elem((i-1)*NNE+j,2)=elem(i,j);
        nw_elem((i-1)*NNE+j,3)=elem(i,1);
        
        nw_mat((i-1)*NNE+j,1)=mat(i);
        
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

function [nw_elem,patch_con,patch_el,nw_material]=split2(x_a,elem,material)

    [elements,NNE]=size(elem);
    [~,sp]=size(x_a);
        
    nw_elem=zeros(elements*sp,NNE-1);
    nw_material=zeros(elements*sp,1);
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
        
        nw_material((i-1)*sp+1,1)=material(i);
        nw_material((i-1)*sp+2,1)=material(i);
        
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

function [elem_0]=corner_nds(elem_1)
    elem_0(:,1)=elem_1(:,1);
    elem_0(:,2)=elem_1(:,3);
    elem_0(:,3)=elem_1(:,5);
    elem_0(:,4)=elem_1(:,7);
    
    elem_0(:,5)=elem_1(:,2);
    elem_0(:,6)=elem_1(:,4);
    elem_0(:,7)=elem_1(:,6);
    elem_0(:,8)=elem_1(:,8);
end

function [elem_0]=corner_nds6(elem_1)
    elem_0(:,1)=elem_1(:,1);
    elem_0(:,2)=elem_1(:,3);
    elem_0(:,3)=elem_1(:,5);
    
    elem_0(:,4)=elem_1(:,2);
    elem_0(:,5)=elem_1(:,4);
    elem_0(:,6)=elem_1(:,6);
end

function [xg,Area,patch_el,patch_con,nw_mat]=quad4xg(x_a,elem_0,mat)


    [~,sp]=size(x_a);
    
    [elements,NNE]=size(elem_0);
    
    [x0]=q_g_center(x_a,elem_0);
    
    xg=zeros(elements*4,sp);
    Area=zeros(elements*4,1);
    nw_mat=zeros(elements*4,1);
    patch_el=zeros(elements*4,1);
    patch_con=zeros(elements,4);
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
            %elem(k,:)=elem_0(e,:);
            Area(k)=area_/NNE;
            nw_mat(k)=mat(e);
            
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

function [xg,Area,patch_el,patch_con,nw_mat]=tri3xg(x_a,elem_0,mat)


    [~,sp]=size(x_a);
    
    [elements,NNE]=size(elem_0);
    
    [~,area_]=g_center(x_a,elem_0,0);
    
    xg=zeros(elements*3,sp);
    Area=zeros(elements*3,1);
    nw_mat=zeros(elements*3,1);
    patch_el=zeros(elements*3,1);
    patch_con=zeros(elements,3);
    
    point=[0.2 0.2; 0.2 0.6; 0.6 0.2];
    
    k=0;
    for e=1:elements
        xn=zeros(NNE,1);
        yn=zeros(NNE,1);

        for i=1:NNE
            nd=elem_0(e,i);
            xn(i)=x_a(nd,1);
            yn(i)=x_a(nd,2);
        end
        
        for i=1:NNE
            
            xp=[xn(1)+(xn(2)-xn(1))*point(i,1)+(xn(3)-xn(1))*point(i,2),...
                yn(1)+(yn(2)-yn(1))*point(i,1)+(yn(3)-yn(1))*point(i,2)];
            
            k=k+1;
            
            for j=1:sp
                xg(k,j)=xp(j);
            end
            %elem(k,:)=elem_0(e,:);
            Area(k)=area_(e)/NNE;
            nw_mat(k)=mat(e);
            
        end 
        for j=1:3
            patch_con(e,j)=(e-1)*NNE+j;
        end
    end
    
    for j=1:elements
        for k=1:NNE
            patch_el(patch_con(j,k))=j;
        end
    end
end

function [xg,Area,nw_elem,patch_con,patch_el,nw_mat]=...
                    split3xg(x_a,elem,mat)
                
    [elements,NNE]=size(elem);
    [~,sp]=size(x_a);
        
    nw_elem=zeros(elements*sp,NNE-1);
    material=zeros(elements*sp,1);

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
        
        material((i-1)*sp+1,1)=mat(i);
        material((i-1)*sp+2,1)=mat(i);
    end
    
    [~,area_]=g_center(x_a,nw_elem,0);
    
    xg=zeros(elements*3*2,sp);
    Area=zeros(elements*3*2,1);
    nw_mat=zeros(elements*3*2,1);
    patch_el=zeros(elements*3*2,1);
    patch_con=zeros(elements,3*2);
    
    point=[0.2 0.2; 0.2 0.6; 0.6 0.2];
    
    k=0;
    for e=1:elements*2
        xn=zeros(3,1);
        yn=zeros(3,1);

        for i=1:3
            nd=nw_elem(e,i);
            xn(i)=x_a(nd,1);
            yn(i)=x_a(nd,2);
        end
        
        for i=1:3
            
            xp=[xn(1)+(xn(2)-xn(1))*point(i,1)+(xn(3)-xn(1))*point(i,2),...
                yn(1)+(yn(2)-yn(1))*point(i,1)+(yn(3)-yn(1))*point(i,2)];
            
            k=k+1;
            
            for j=1:sp
                xg(k,j)=xp(j);
            end
            %elem(k,:)=elem_0(e,:);
            Area(k)=area_(e)/3;
            nw_mat(k)=material(e);
            
        end 
    end
    
    for e=1:elements
        for j=1:3*2
            patch_con(e,j)=(e-1)*6+j;
        end
    end
    
    for j=1:elements
        for k=1:3*sp
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

function [x,elem,NNE,material,NODE_LIST]=read_dat(str1,str2,DIM)
%

    close all
    
    % File: Geo_DynCLM
    %   Read geometry in *.dat and save in the geometry variables

    % Date:
    %   Version 3.0   25.11.2019

    %Input file name
    if ismac || isunix  % Code to run on Mac or Unix plaform 
        filename_i=strcat(str2,'/',str1,'.dat');
    elseif ispc         % Code to run on Windows platform
        filename_i=strcat(str2,'\',str1,'.dat');
    else
        disp('Platform not supported')
        stop
    end 

    f_i = fopen(filename_i, 'rt');
    formato = '%s %s %s %s %s %s %s %s %s %s %s'; % formato de cada línea 

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
    
    [l,~]=size(data);
    
    t=0;
    while (t<l)
        t=t+1;
        s1=data{t,1};

        switch s1
            case 'NNODE='
                nodes=str2double(data{t,2});
                continue
            case 'NCONE='
                NNE=str2double(data{t,2});
                continue
            case 'NDIME='
                if DIM>str2double(data{t,2})
                    fprintf('Error, DIM mismatches with spatial dimension of the mesh!!\n')
                    stop
                end
                x=zeros(nodes,DIM);
                continue
            case 'NELEM='
                elements=str2double(data{t,2});
                elem=zeros(elements,NNE);
                material=zeros(elements,1);
                continue
            case 'NMATS='
                mats=str2double(data{t,2});
                continue
            case 'CONNECTIVITY'
                break
            otherwise
                continue
        end
    end
    
    % CHECK
    if t>l
        fprintf('Error, bad reading of geometry!!\n')
        stop
    end
    
    % CONNECTIVITY
    i=0;
    for e=1:elements
        t=t+1;
        if str2double(data{t,2})~=0
            i=i+1;
            material(i)=str2double(data{t,2});
            for j=1:NNE
                elem(i,j)=str2double(data{t,2+j});
            end
        end
    end
    elem=elem(1:i,:);
    
    % CHECK
    t=t+2;
    if strcmp(data{t,1},'NODES')~=1
        fprintf('Error, bad reading of geometry!!\n')
        stop
    end
    
    % NODES
    for i=1:nodes
        t=t+1;
        for j=1:DIM
            x(i,j)=str2double(data{t,1+j});
        end
    end

    t=t+4;
    fin=0;
    if t>l
        fin=1;
    elseif strcmp(data{t,1},'BOUNDARY_CONDITIONS')==0 && ...
            strcmp(data{t,1},'LOADS')==0
        fprintf('Error, bad reading of geometry!!\n')
        stop
    end
    
    if fin==0
        bcs=0;
        pls=0;
        lls=0;
        vls=0;
        abcs=0;
        rbs=0;
        rbs_nds=[];
        BC={};
        PL={};
        ABC={};
        LL={};
        VL={};
        RB={};
    end
    
    while fin==0
        t=t+1;
        s1=data{t,1};
        switch s1
            case 'BC_SET'
                num=str2double(data{t,3});
                if num~=0
                    bcs=bcs+1;
                    if bcs==str2double(data{t,2})
                        list=zeros(num,1);
                        for i=1:num
                            t=t+1;
                            list(i)=str2double(data{t,1});
                        end
                        BC{bcs}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of BC!!\n')
                        stop
                    end
                end
                continue
            case 'RIGID_BODY'
                num=str2double(data{t,3});
                if num~=0
                    rbs=rbs+1;
                    if rbs==str2double(data{t,2})
                        list=zeros(num,3);
                        for i=1:num
                            t=t+1;
                            list(i,1)=x(str2double(data{t,1}),1);
                            list(i,2)=x(str2double(data{t,1}),2);
                            list(i,3)=x(str2double(data{t,2}),1);
                            list(i,4)=x(str2double(data{t,2}),2);
                        end
                        rbs_nds=union(rbs_nds,[str2double(data{t,1}) str2double(data{t,2})]);
                        RB{rbs}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of LOADs!!\n')
                        stop
                    end
                end
                continue
            case 'ABSORBING_BC'
                num=str2double(data{t,3});
                if num~=0
                    abcs=abcs+1;
                    if abcs==str2double(data{t,2})
                        list=zeros(num,3);
                        for i=1:num
                            t=t+1;
                            list(i,1)=str2double(data{t,1});
                            list(i,2)=str2double(data{t,2});
                            if NNE==8 || NNE==6
                              list(i,3)=str2double(data{t,3});  
                            end
                        end
                        ABC{abcs}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of LOADs!!\n')
                        stop
                    end
                end
                continue
            case 'POINT_LOAD'
                num=str2double(data{t,3});
                if num~=0
                    pls=pls+1;
                    if pls==str2double(data{t,2})
                        list=zeros(num,1);
                        for i=1:num
                            t=t+1;
                            list(i)=str2double(data{t,1});
                        end
                        PL{pls}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of LOADs!!\n')
                        stop
                    end
                end
                continue
            case 'LINE_LOAD'
                num=str2double(data{t,3});
                if num~=0
                    lls=lls+1;
                    if lls==str2double(data{t,2})
                        list=zeros(num,3);
                        for i=1:num
                            t=t+1;
                            list(i,1)=str2double(data{t,1});
                            list(i,2)=str2double(data{t,2});
                            if NNE==8 || NNE==6
                              list(i,3)=str2double(data{t,3});  
                            end
                        end
                        LL{lls}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of LOADs!!\n')
                        stop
                    end
                end
                continue
            case 'VOLUME_LOAD'
                num=str2double(data{t,3});
                if num~=0
                    vls=vls+1;
                    if vls==str2double(data{t,2})
                        list=zeros(num,1);
                        for i=1:num
                            t=t+1;
                            list(i)=str2double(data{t,1});
                        end
                        VL{vls}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of LOADs!!\n')
                        stop
                    end 
                end
                continue
            case 'END_BC'
                t=t+2;
                if t>l
                    fin=1;
                elseif strcmp(data{t,1},'LOADS')==0
                    fprintf('Error, bad reading of geometry!!\n')
                    stop
                end
                continue
            case 'END_LDS'
                fin=1;
                continue
            otherwise
                fprintf('Error, bad reading of geometry!!\n')
                stop
        end
    end
    
    if isempty(rbs_nds)==0
       [a,b]=size(elem);
       for i=length(rbs_nds):-1:1
           x=[x(1:rbs_nds(i)-1,:);x(rbs_nds(i)+1:end,:)];
           for j=1:a
               for k=1:b
                   if elem(j,k)>rbs_nds(i)
                       elem(j,k)=elem(j,k)-1;
                   end
               end
           end
           %BCS
           for j=1:bcs
               list=BC{j};
               for k=1:length(list)
                   if list(k)>rbs_nds(i)
                        list(k)=list(k)-1;
                   end
               end 
               BC{j}=list;
           end
           %PLS
           for j=1:pls
               list=PL{j};
               for k=1:length(list)
                   if list(k)>rbs_nds(i)
                        list(k)=list(k)-1;
                   end
               end
               PL{j}=list;
           end
           %VLS
           for j=1:vls
               list=VL{j};
               for k=1:length(list)
                   if list(k)>rbs_nds(i)
                        list(k)=list(k)-1;
                   end
               end  
               VL{j}=list;
           end
           %LLS
           for j=1:lls
               list=LL{j};
               [a1,b1]=size(list);
               for k=1:a1
                   for l=1:b1
                        if list(k,l)>rbs_nds(i)
                            list(k,l)=list(k,l)-1;
                        end
                   end
               end
               LL{j}=list;
           end
       end
    end

    NODE_LIST=struct('bcs',bcs,'pls',pls,'lls',lls,'vls',vls,'rbs',rbs,...
        'abcs',abcs,'BC',{BC},'PL',{PL},'LL',{LL},'VL',{VL},'RB',{RB},...
        'ABC',{ABC});
end

