function vtk_driver(str,str2,steps,h,rt,folder)

    root=strcat(folder,'/',str,'.mat');

    load(root,'GEOMETRY','SOLVER');

    [~,NNE]=size(GEOMETRY.elem);
    
    if strcmp(steps,'FULL')
        steps=SOLVER.dim;
    end

    if NNE==4 || NNE==8
        if SOLVER.UW==1
            vtk_4DOF_quad(str,str2,steps,h,rt);
        elseif SOLVER.UW==2
            vtk_3DOF_quad(str,str2,steps,h,rt);
        elseif SOLVER.UW==0
            vtk_2DOF_quad(root,str2,steps,h,rt);
        end
    else
        if SOLVER.UW==1
            vtk_4DOF(str,str2,steps,h,rt);
        elseif SOLVER.UW==2
            vtk_3DOF(str,str2,steps,h,rt);
        elseif SOLVER.UW==0
            vtk_2DOF(str,str2,steps,h,rt);
        end
    end


end


function vtk_2DOF(str,str2,steps,sc)

load(str,'-mat');
[elements,NNE]=size(elem);
[nodes,sp]=size(x_a);
df=sp;


dd=round(ste_p/steps);

x_a=x_0;

%[~,ste_p]=size(Ps);

for cont=1:dd:ste_p

    %Output file name
    filename=(['VTK/' str2 '_' num2str(cont) '.vtk']);
   

    %nr_of_elements=numel(x);
    fid = fopen(filename, 'w'); 

    %ASCII file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'vtk output\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, ['POINTS ' num2str(nodes) ' float\n']);

    for i=1:nodes
        fprintf(fid,[num2str(x_a(i,1)+sc*d(i*df-1,cont)) ' ' num2str(x_a(i,2)+sc*d(i*df,cont)) ' 0\n']);
    end
    %grid	
    fprintf(fid, ['CELLS ' num2str(elements) ' ' num2str(elements*(NNE+1)) '\n']);
    for i=1:elements
        fprintf(fid,[num2str(NNE) ' ' num2str(elem(i,1)-1) ' ' num2str(elem(i,2)-1) ' ' num2str(elem(i,3)-1) '\n']);
    end
    fprintf(fid, ['CELL_TYPES ' num2str(elements) '\n']);
    for i=1:elements
        fprintf(fid,' 5 \n');
    end
    fprintf(fid, ['POINT_DATA ' num2str(nodes) '\n']);
    %d
    fprintf(fid, 'VECTORS d float \n');
    for i=1:nodes
        fprintf(fid,[num2str(d(i*df-1,cont)) ' ' num2str(d(i*df,cont)) ' 0\n']);
    end
    %V
    fprintf(fid, 'VECTORS v float \n');
    for i=1:nodes
        fprintf(fid,[num2str(v(i*df-1,cont)) ' ' num2str(v(i*df,cont)) ' 0\n']);
    end
    %Ep
    fprintf(fid, 'SCALARS gamma float \n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:nodes
        fprintf(fid,[num2str(Gamma_nds(i,cont)) '\n']);
    end
    %Aw
%     fprintf(fid, 'VECTORS a_w float \n');
%     for i=1:nodes
%         fprintf(fid,[num2str(v(i*df-1,cont)) ' ' num2str(a(i*df,cont)) ' 0\n']);
%     end
    
    fprintf(fid, ['CELL_DATA ' num2str(elements) '\n']);
    %Pressure
    fprintf(fid,'SCALARS Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ps(i,cont)) '\n']);
    end
    %Sy
    fprintf(fid,'SCALARS Yield_str float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Sy_tot(i,cont)) '\n']);
    end
     %Plastic strain
     fprintf(fid,'SCALARS Ep float\n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:elements
         fprintf(fid,[num2str(Gamma_tot(i,cont)) '\n']);
     end
    
    %E11
    fprintf(fid,'SCALARS E11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-3,cont)) '\n']);
    end
    %E22
    fprintf(fid,'SCALARS E22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-2,cont)) '\n']);
    end
    %E33
    fprintf(fid,'SCALARS E33 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-1,cont)) '\n']);
    end
    %G12
    fprintf(fid,'SCALARS G12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4,cont)) '\n']);
    end
    
    %E11_p
    fprintf(fid,'SCALARS E11_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-3,cont)) '\n']);
    end
    %E22_p
    fprintf(fid,'SCALARS E22_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-2,cont)) '\n']);
    end
    %E33_p
    fprintf(fid,'SCALARS E33_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-1,cont)) '\n']);
    end
    %G12_p
    fprintf(fid,'SCALARS G12_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4,cont)) '\n']);
    end
    
    %S11
    fprintf(fid,'SCALARS S11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-3,cont)) '\n']);
    end
    %S22
    fprintf(fid,'SCALARS S22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-2,cont)) '\n']);
    end
    %S33
    fprintf(fid,'SCALARS S33 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-1,cont)) '\n']);
    end
    %T12
    fprintf(fid,'SCALARS T12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4,cont)) '\n']);
    end

    fclose(fid);
end

end

function vtk_2DOF_quad(str,str2,steps,sc,rt)

    %clear
    load(str,'-mat','GLOBAL','GEOMETRY');
    if strcmp(GEOMETRY.ELEMENT,'Q8P4-4') || strcmp(GEOMETRY.ELEMENT,'Q4-4') ...
        || strcmp(GEOMETRY.ELEMENT,'Q8-4')
        it=4;
    else
        it=1;
    end
    [elements,NNE]=size(GEOMETRY.elem);
    
    x_a=GEOMETRY.x_0;
    [nodes,sp]=size(x_a);
    df=sp;
    %[ste_p,~]=size(tp);
    
    if NNE==4
        elem2=GEOMETRY.elem;
        nodes_p = linspace(1,GEOMETRY.nodes)';
    else
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
    end
    Ss=GLOBAL.Sigma;
    Es=GLOBAL.Es;
    Es_p=GLOBAL.Es_p;
    d=GLOBAL.d;
    a=GLOBAL.a;
    v=GLOBAL.v;

    for cont=1:rt:steps

        %Output file name
        filename=(['VTK/' str2 '_' num2str(cont) '.vtk']);

        if cont==41
            cont;
        end
        %nr_of_elements=numel(x);
        fid = fopen(filename, 'w'); 

        %ASCII file header
        fprintf(fid, '# vtk DataFile Version 3.0\n');
        fprintf(fid, 'vtk output\n');
        fprintf(fid, 'ASCII\n');
        fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
        fprintf(fid, ['POINTS ' num2str(length(nodes_p)) ' float\n']);

        for j=1:length(nodes_p)
            i=nodes_p(j);
            fprintf(fid,[num2str(x_a(i,1)+sc*d(i*df-1,cont),'%10.5f\t') ' ' num2str(x_a(i,2)+sc*d(i*df,cont),'%10.5f\t') ' 0\n']);
        end
        %grid	
        fprintf(fid, ['CELLS ' num2str(elements) ' ' num2str(elements*(NNE+1)) '\n']);
        for i=1:elements
            fprintf(fid,[num2str(NNE) ' ' num2str(elem2(i,1)-1) ' ' num2str(elem2(i,2)-1) ' ' num2str(elem2(i,3)-1) ' ' num2str(elem2(i,4)-1) '\n']);
        end
        fprintf(fid, ['CELL_TYPES ' num2str(elements) '\n']);
        for i=1:elements
            fprintf(fid,' 9 \n');
        end
        fprintf(fid, ['POINT_DATA ' num2str(length(nodes_p)) '\n']);
        %d
        fprintf(fid, 'VECTORS d float \n');
        for j=1:length(nodes_p)
            i=nodes_p(j);
            fprintf(fid,[num2str(d(i*df-1,cont),'%20.10f\t') ' ' num2str(d(i*df,cont),'%20.10f\t') ' 0\n']);
        end
        %V
        fprintf(fid, 'VECTORS v float \n');
        for j=1:length(nodes_p)
            i=nodes_p(j);
            fprintf(fid,[num2str(v(i*df-1,cont),'%20.10f\t') ' ' num2str(v(i*df,cont),'%20.10f\t') ' 0\n']);
        end
        %A
        fprintf(fid, 'VECTORS a float \n');
        for j=1:length(nodes_p)
            i=nodes_p(j);
            fprintf(fid,[num2str(a(i*df-1,cont),'%20.10f\t') ' ' num2str(a(i*df,cont),'%20.10f\t') ' 0\n']);
        end
        %Aw
    %     fprintf(fid, 'VECTORS a_w float \n');
    %     for i=1:nodes
    %         fprintf(fid,[num2str(v(i*df-1,cont)) ' ' num2str(a(i*df,cont)) ' 0\n']);
    %     end

        fprintf(fid, ['CELL_DATA ' num2str(elements) '\n']);
        %Pressure
        fprintf(fid,'SCALARS Pressure float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(GLOBAL.Ps(i*it,cont),'%20.10f\t') '\n']);
        end
        %Q
        fprintf(fid,'SCALARS Q float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(GLOBAL.Qs(i*it,cont),'%20.10f\t') '\n']);
        end
        %Sy
        fprintf(fid,'SCALARS Yield_str float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(GLOBAL.Sy(i*it,cont),'%20.10f\t') '\n']);
        end
         %Plastic strain
         fprintf(fid,'SCALARS Ep float\n');
         fprintf(fid, 'LOOKUP_TABLE default\n');
         for i=1:elements
             fprintf(fid,[num2str(GLOBAL.gamma(i*it,cont),'%20.10f\t') '\n']);
         end

        %E11
        fprintf(fid,'SCALARS E11 float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Es(i*it*4-3,cont),'%20.10f\t') '\n']);
        end
        %E22
        fprintf(fid,'SCALARS E22 float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Es(i*it*4-2,cont),'%20.10f\t') '\n']);
        end
        %E33
        fprintf(fid,'SCALARS E33 float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Es(i*it*4-1,cont),'%20.10f\t') '\n']);
        end
        %G12
        fprintf(fid,'SCALARS G12 float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Es(i*it*4,cont),'%20.10f\t') '\n']);
        end

        %E11_p
        fprintf(fid,'SCALARS E11_p float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Es_p(i*it*4-3,cont),'%20.10f\t') '\n']);
        end
        %E22_p
        fprintf(fid,'SCALARS E22_p float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Es_p(i*it*4-2,cont),'%20.10f\t') '\n']);
        end
        %E33_p
        fprintf(fid,'SCALARS E33_p float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Es_p(i*it*4-1,cont),'%20.10f\t') '\n']);
        end
        %G12_p
        fprintf(fid,'SCALARS G12_p float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Es_p(i*it*4,cont),'%20.10f\t') '\n']);
        end

        %S11
        fprintf(fid,'SCALARS S11 float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Ss(i*it*4-3,cont),'%20.10f\t') '\n']);
        end
        %S22
        fprintf(fid,'SCALARS S22 float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Ss(i*it*4-2,cont),'%20.10f\t') '\n']);
        end
        %S33
        fprintf(fid,'SCALARS S33 float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Ss(i*it*4-1,cont),'%20.10f\t') '\n']);
        end
        %T12
        fprintf(fid,'SCALARS T12 float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Ss(i*it*4,cont),'%20.10f\t') '\n']);
        end

        fclose(fid);
    end

end
 
function vtk_4DOF(str,str2,steps,sc,dd)

%clear
load(str,'-mat','GLOBAL','GEOMETRY');
[elements,NNE]=size(GEOMETRY.elem);
x_a=GEOMETRY.x_0;
[nodes,sp]=size(x_a);
df=2*sp;
%[ste_p,~]=size(tp);

elem=GEOMETRY.elem;
Ss=GLOBAL.Sigma;
Es=GLOBAL.Es;
Pw=GLOBAL.pw;
Ps=GLOBAL.Ps;
Qs=GLOBAL.Qs;
Es_p=GLOBAL.Es_p;
Sy_tot=GLOBAL.Sy;
gamma_nds=GLOBAL.gamma_nds;
gamma=GLOBAL.gamma;
d=GLOBAL.d;
a=GLOBAL.a;
v=GLOBAL.v;

for cont=1:dd:steps

    %Output file name
    filename=(['VTK/' str2 '_' num2str(cont) '.vtk']);
   

    %nr_of_elements=numel(x);
    fid = fopen(filename, 'w'); 

    %ASCII file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'vtk output\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, ['POINTS ' num2str(nodes) ' float\n']);

    for i=1:nodes
        fprintf(fid,[num2str(x_a(i,1)+sc*d(i*df-3,cont)) ' ' num2str(x_a(i,2)+sc*d(i*df-2,cont)) ' 0\n']);
    end
    %grid	
    fprintf(fid, ['CELLS ' num2str(elements) ' ' num2str(elements*(NNE+1)) '\n']);
    for i=1:elements
        fprintf(fid,[num2str(NNE) ' ' num2str(elem(i,1)-1) ' ' num2str(elem(i,2)-1) ' ' num2str(elem(i,3)-1) '\n']);
    end
    fprintf(fid, ['CELL_TYPES ' num2str(elements) '\n']);
    for i=1:elements
        fprintf(fid,' 5 \n');
    end
    fprintf(fid, ['POINT_DATA ' num2str(nodes) '\n']);
    %dw
    fprintf(fid, 'VECTORS d_w float \n');
    for i=1:nodes
        fprintf(fid,[num2str(d(i*df-1,cont)) ' ' num2str(d(i*df,cont)) ' 0\n']);
    end
    %du
    fprintf(fid, 'VECTORS d_u float \n');
    for i=1:nodes
        fprintf(fid,[num2str(d(i*df-3,cont)) ' ' num2str(d(i*df-2,cont)) ' 0\n']);
    end
    %Vw
    fprintf(fid, 'VECTORS v_w float \n');
    for i=1:nodes
        fprintf(fid,[num2str(v(i*df-1,cont)) ' ' num2str(v(i*df,cont)) ' 0\n']);
    end
    %Vu
    fprintf(fid, 'VECTORS v_u float \n');
    for i=1:nodes
        fprintf(fid,[num2str(v(i*df-3,cont)) ' ' num2str(v(i*df-2,cont)) ' 0\n']);
    end
    %Aw
     fprintf(fid, 'VECTORS a_w float \n');
     for i=1:nodes
         fprintf(fid,[num2str(a(i*df-1,cont)) ' ' num2str(a(i*df,cont)) ' 0\n']);
     end
     %Au
     fprintf(fid, 'VECTORS a_u float \n');
     for i=1:nodes
         fprintf(fid,[num2str(a(i*df-3,cont)) ' ' num2str(a(i*df-2,cont)) ' 0\n']);
     end
     %Gamma
     fprintf(fid, 'SCALARS Ep_nds float \n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:nodes
         fprintf(fid,[num2str(gamma_nds(i,cont)) '\n']);
     end

    
    fprintf(fid, ['CELL_DATA ' num2str(elements) '\n']);
    %Pressure
    fprintf(fid,'SCALARS Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ps(i,cont)) '\n']);
    end
    %Q
    fprintf(fid,'SCALARS Q float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Qs(i,cont)) '\n']);
    end
    
    %Sy
    fprintf(fid,'SCALARS Yield_str float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Sy_tot(i,cont)) '\n']);
    end
     %Plastic strain
     fprintf(fid,'SCALARS Ep float\n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:elements
         fprintf(fid,[num2str(gamma(i,cont)) '\n']);
     end
    %Pore Pressure
    fprintf(fid,'SCALARS Pore_Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Pw(i,cont)) '\n']);
    end
    
    %E11
    fprintf(fid,'SCALARS E11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-3,cont)) '\n']);
    end
    %E22
    fprintf(fid,'SCALARS E22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-2,cont)) '\n']);
    end
    %E22
    fprintf(fid,'SCALARS E33 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-1,cont)) '\n']);
    end
    %G12
    fprintf(fid,'SCALARS G12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4,cont)) '\n']);
    end
%     %EW11
%     fprintf(fid,'SCALARS EW11 float\n');
%     fprintf(fid, 'LOOKUP_TABLE default\n');
%     for i=1:elements
%         fprintf(fid,[num2str(Es_w(i*5-1,cont)) '\n']);
%     end
%     %EW22
%     fprintf(fid,'SCALARS EW22 float\n');
%     fprintf(fid, 'LOOKUP_TABLE default\n');
%     for i=1:elements
%         fprintf(fid,[num2str(Es_w(i*4,cont)) '\n']);
%     end

    %E11_p
    fprintf(fid,'SCALARS E11_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-3,cont)) '\n']);
    end
    %E22_p
    fprintf(fid,'SCALARS E22_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-2,cont)) '\n']);
    end
    %E33_p
    fprintf(fid,'SCALARS E33_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-1,cont)) '\n']);
    end
    %G12_p
    fprintf(fid,'SCALARS G12_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4,cont)) '\n']);
    end
    
    %S11
    fprintf(fid,'SCALARS S11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-3,cont)) '\n']);
    end
    %S22
    fprintf(fid,'SCALARS S22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-2,cont)) '\n']);
    end
    %S33
    fprintf(fid,'SCALARS S33 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-1,cont)) '\n']);
    end
    %T12
    fprintf(fid,'SCALARS T12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4,cont)) '\n']);
    end


    fclose(fid);
end

end

function vtk_4DOF_quad(str,str2,steps,sc,dd)

%clear
load(str,'-mat','GLOBAL','GEOMETRY');
if strcmp(GEOMETRY.ELEMENT,'Q8P4-4') || strcmp(GEOMETRY.ELEMENT,'Q4-4') ...
        || strcmp(GEOMETRY.ELEMENT,'Q8-4')
    it=4;
else
    it=1;
end
[elements,NNE]=size(GEOMETRY.elem);
x_a=GEOMETRY.x_0;
[nodes,sp]=size(x_a);
df=2*sp;
%[ste_p,~]=size(tp);

    if NNE==4
        elem=GEOMETRY.elem;
        elem2=GEOMETRY.elem;
        nodes_p = linspace(1,GEOMETRY.nodes)';
    else
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
    end
Ss=GLOBAL.Sigma;
Es=GLOBAL.Es;
Pw=GLOBAL.pw;
Ps=GLOBAL.Ps;
Es_p=GLOBAL.Es_p;
Sy_tot=GLOBAL.Sy;
gamma_nds=GLOBAL.gamma_nds;
gamma=GLOBAL.gamma;
d=GLOBAL.d;
a=GLOBAL.a;
v=GLOBAL.v;

for cont=1:dd:steps

    %Output file name
    filename=(['VTK/' str2 '_' num2str(cont) '.vtk']);
   

    %nr_of_elements=numel(x);
    fid = fopen(filename, 'w'); 

    %ASCII file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'vtk output\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, ['POINTS ' num2str(length(nodes_p)) ' float\n']);

    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(x_a(i,1)+sc*d(i*df-3,cont)) ' ' num2str(x_a(i,2)+sc*d(i*df-2,cont)) ' 0\n']);
    end
    %grid	
    fprintf(fid, ['CELLS ' num2str(elements) ' ' num2str(elements*(NNE+1)) '\n']);
    for i=1:elements
        fprintf(fid,[num2str(NNE) ' ' num2str(elem2(i,1)-1) ' ' num2str(elem2(i,2)-1) ' ' num2str(elem2(i,3)-1) ' ' num2str(elem2(i,4)-1) '\n']);
    end
    fprintf(fid, ['CELL_TYPES ' num2str(elements) '\n']);
    for i=1:elements
        fprintf(fid,' 9 \n');
    end
    fprintf(fid, ['POINT_DATA ' num2str(length(nodes_p)) '\n']);
    %dw
    fprintf(fid, 'VECTORS d_w float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(d(i*df-1,cont)) ' ' num2str(d(i*df,cont)) ' 0\n']);
    end
    %du
    fprintf(fid, 'VECTORS d_u float \n');
     for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(d(i*df-3,cont)) ' ' num2str(d(i*df-2,cont)) ' 0\n']);
    end
    %Vw
    fprintf(fid, 'VECTORS v_w float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(v(i*df-1,cont)) ' ' num2str(v(i*df,cont)) ' 0\n']);
    end
    %Vu
    fprintf(fid, 'VECTORS v_u float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(v(i*df-3,cont)) ' ' num2str(v(i*df-2,cont)) ' 0\n']);
    end
    %Aw
     fprintf(fid, 'VECTORS a_w float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
         fprintf(fid,[num2str(a(i*df-1,cont)) ' ' num2str(a(i*df,cont)) ' 0\n']);
     end
     %Au
     fprintf(fid, 'VECTORS a_u float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
         fprintf(fid,[num2str(a(i*df-3,cont)) ' ' num2str(a(i*df-2,cont)) ' 0\n']);
     end
     %Gamma
     fprintf(fid, 'SCALARS Ep_nds float \n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
         fprintf(fid,[num2str(gamma_nds(i,cont)) '\n']);
     end

    
    fprintf(fid, ['CELL_DATA ' num2str(elements) '\n']);
    %Pressure
    fprintf(fid,'SCALARS Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ps(i*it,cont)) '\n']);
    end
    %Sy
    fprintf(fid,'SCALARS Yield_str float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Sy_tot(i*it,cont)) '\n']);
    end
     %Plastic strain
     fprintf(fid,'SCALARS Ep float\n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:elements
         fprintf(fid,[num2str(gamma(i*it,cont)) '\n']);
     end
    %Pore Pressure
    fprintf(fid,'SCALARS Pore_Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Pw(i*it,cont)) '\n']);
    end
    
    %E11
    fprintf(fid,'SCALARS E11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*it*4-3,cont)) '\n']);
    end
    %E22
    fprintf(fid,'SCALARS E22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*it*4-2,cont)) '\n']);
    end
    %G12
    fprintf(fid,'SCALARS G12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*it*4,cont)) '\n']);
    end
%     %EW11
%     fprintf(fid,'SCALARS EW11 float\n');
%     fprintf(fid, 'LOOKUP_TABLE default\n');
%     for i=1:elements
%         fprintf(fid,[num2str(Es_w(i*5-1,cont)) '\n']);
%     end
%     %EW22
%     fprintf(fid,'SCALARS EW22 float\n');
%     fprintf(fid, 'LOOKUP_TABLE default\n');
%     for i=1:elements
%         fprintf(fid,[num2str(Es_w(i*4,cont)) '\n']);
%     end

    fclose(fid);
end

end
 
function vtk_3DOF(str,str2,steps,sc,dd)

%clear
load(str,'-mat','GLOBAL','GEOMETRY');
[elements,NNE]=size(GEOMETRY.elem);
x_a=GEOMETRY.x_0;
[nodes,sp]=size(x_a);
df=sp+1;
%[ste_p,~]=size(tp);

elem=GEOMETRY.elem;
Ss=GLOBAL.Sigma;
Es=GLOBAL.Es;
Pw=GLOBAL.pw;
dPw=GLOBAL.dpw;
Ps=GLOBAL.Ps;
Qs=GLOBAL.Qs;
Es_p=GLOBAL.Es_p;
Sy_tot=GLOBAL.Sy;
gamma_nds=GLOBAL.gamma_nds;
gamma=GLOBAL.gamma;
d=GLOBAL.d;
a=GLOBAL.a;
v=GLOBAL.v;

for cont=1:dd:steps

    %Output file name
    filename=(['VTK/' str2 '_' num2str(cont) '.vtk']);
   

    %nr_of_elements=numel(x);
    fid = fopen(filename, 'w'); 

    %ASCII file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'vtk output\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, ['POINTS ' num2str(nodes) ' float\n']);

    for i=1:nodes
        fprintf(fid,[num2str(x_a(i,1)+sc*d(i*df-2,cont)) ' ' num2str(x_a(i,2)+sc*d(i*df-1,cont)) ' 0\n']);
    end
    %grid	
    fprintf(fid, ['CELLS ' num2str(elements) ' ' num2str(elements*(NNE+1)) '\n']);
    for i=1:elements
        fprintf(fid,[num2str(NNE) ' ' num2str(elem(i,1)-1) ' ' num2str(elem(i,2)-1) ' ' num2str(elem(i,3)-1) '\n']);
    end
    fprintf(fid, ['CELL_TYPES ' num2str(elements) '\n']);
    for i=1:elements
        fprintf(fid,' 5 \n');
    end
    fprintf(fid, ['POINT_DATA ' num2str(nodes) '\n']);
    %pw
    fprintf(fid, 'SCALARS pw_nds float \n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:nodes
        fprintf(fid,[num2str(d(i*df,cont)) '\n']);
    end
    %vpw
    fprintf(fid, 'SCALARS vpw_nds float \n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:nodes
        fprintf(fid,[num2str(v(i*df,cont)) '\n']);
    end
      %apw
    fprintf(fid, 'SCALARS apw_nds float \n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:nodes
         fprintf(fid,[num2str(a(i*df,cont)) '\n']);
     end
    
    %du
    fprintf(fid, 'VECTORS d_u float \n');
    for i=1:nodes
        fprintf(fid,[num2str(d(i*df-2,cont)) ' ' num2str(d(i*df-1,cont)) ' 0\n']);
    end
    %Vu
    fprintf(fid, 'VECTORS v_u float \n');
    for i=1:nodes
        fprintf(fid,[num2str(v(i*df-2,cont)) ' ' num2str(v(i*df-1,cont)) ' 0\n']);
    end
     %Au
     fprintf(fid, 'VECTORS a_u float \n');
     for i=1:nodes
         fprintf(fid,[num2str(a(i*df-2,cont)) ' ' num2str(a(i*df-1,cont)) ' 0\n']);
     end
     
     %Gamma
     fprintf(fid, 'SCALARS Ep_nds float \n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:nodes
         fprintf(fid,[num2str(gamma_nds(i,cont)) '\n']);
     end

    
    fprintf(fid, ['CELL_DATA ' num2str(elements) '\n']);
    %Pressure
    fprintf(fid,'SCALARS Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ps(i,cont)) '\n']);
    end
    %Q
    fprintf(fid,'SCALARS Q float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Qs(i,cont)) '\n']);
    end
    
    %Sy
    fprintf(fid,'SCALARS Yield_str float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Sy_tot(i,cont)) '\n']);
    end
     %Plastic strain
     fprintf(fid,'SCALARS Ep float\n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:elements
         fprintf(fid,[num2str(gamma(i,cont)) '\n']);
     end
    %Pore Pressure
    fprintf(fid,'SCALARS Pore_Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Pw(i,cont)) '\n']);
    end
    
    %Grad Pore Pressure
    fprintf(fid, 'VECTORS Grad_Pore_Pressure float \n');
    for i=1:elements
        fprintf(fid,[num2str(dPw((i-1)*sp+1,cont)) ' ' num2str(dPw((i-1)*sp+2,cont)) ' 0\n']);
    end
    
    %E11
    fprintf(fid,'SCALARS E11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-3,cont)) '\n']);
    end
    %E22
    fprintf(fid,'SCALARS E22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-2,cont)) '\n']);
    end
    %E22
    fprintf(fid,'SCALARS E33 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4-1,cont)) '\n']);
    end
    %G12
    fprintf(fid,'SCALARS G12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*4,cont)) '\n']);
    end
%     %EW11
%     fprintf(fid,'SCALARS EW11 float\n');
%     fprintf(fid, 'LOOKUP_TABLE default\n');
%     for i=1:elements
%         fprintf(fid,[num2str(Es_w(i*5-1,cont)) '\n']);
%     end
%     %EW22
%     fprintf(fid,'SCALARS EW22 float\n');
%     fprintf(fid, 'LOOKUP_TABLE default\n');
%     for i=1:elements
%         fprintf(fid,[num2str(Es_w(i*4,cont)) '\n']);
%     end

    %E11_p
    fprintf(fid,'SCALARS E11_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-3,cont)) '\n']);
    end
    %E22_p
    fprintf(fid,'SCALARS E22_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-2,cont)) '\n']);
    end
    %E33_p
    fprintf(fid,'SCALARS E33_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4-1,cont)) '\n']);
    end
    %G12_p
    fprintf(fid,'SCALARS G12_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*4,cont)) '\n']);
    end
    
    %S11
    fprintf(fid,'SCALARS S11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-3,cont)) '\n']);
    end
    %S22
    fprintf(fid,'SCALARS S22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-2,cont)) '\n']);
    end
    %S33
    fprintf(fid,'SCALARS S33 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4-1,cont)) '\n']);
    end
    %T12
    fprintf(fid,'SCALARS T12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*4,cont)) '\n']);
    end


    fclose(fid);
end

end

function vtk_3DOF_quad(str,str2,steps,sc,dd)

%clear
load(str,'-mat','GLOBAL','GEOMETRY');
if strcmp(GEOMETRY.ELEMENT,'Q8P4-4') || strcmp(GEOMETRY.ELEMENT,'Q4-4') ...
        || strcmp(GEOMETRY.ELEMENT,'Q8-4')
    it=4;
else
    it=1;
end
[elements,NNE]=size(GEOMETRY.elem);
x_a=GEOMETRY.x_0;
[nodes,sp]=size(x_a);
df=sp+1;
%[ste_p,~]=size(tp);

if NNE==4
    elem=GEOMETRY.elem;
    elem2=GEOMETRY.elem;
    nodes_p = linspace(1,GEOMETRY.nodes)';
else
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
end
Ss=GLOBAL.Sigma;
Es=GLOBAL.Es;
Pw=GLOBAL.pw;
dPw=GLOBAL.dpw;
Ps=GLOBAL.Ps;
Es_p=GLOBAL.Es_p;
Sy_tot=GLOBAL.Sy;
gamma_nds=GLOBAL.gamma_nds;
gamma=GLOBAL.gamma;
d=GLOBAL.d;
a=GLOBAL.a;
v=GLOBAL.v;

for cont=1:dd:steps

    %Output file name
    filename=(['VTK/' str2 '_' num2str(cont) '.vtk']);
   
    %nr_of_elements=numel(x);
    fid = fopen(filename, 'w'); 

    %ASCII file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'vtk output\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, ['POINTS ' num2str(length(nodes_p)) ' float\n']);

    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(x_a(i,1)+sc*d(i*df-2,cont)) ' ' num2str(x_a(i,2)+sc*d(i*df-1,cont)) ' 0\n']);
    end
    %grid	
    fprintf(fid, ['CELLS ' num2str(elements) ' ' num2str(elements*(NNE+1)) '\n']);
    for i=1:elements
        fprintf(fid,[num2str(NNE) ' ' num2str(elem2(i,1)-1) ' ' num2str(elem2(i,2)-1) ' ' num2str(elem2(i,3)-1) ' ' num2str(elem2(i,4)-1) '\n']);
    end
    fprintf(fid, ['CELL_TYPES ' num2str(elements) '\n']);
    for i=1:elements
        fprintf(fid,' 9 \n');
    end
    fprintf(fid, ['POINT_DATA ' num2str(length(nodes_p)) '\n']);
        
    %pw
    fprintf(fid, 'SCALARS pw_nds float \n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(d(i*df,cont)) '\n']);
    end
    %vpw
    fprintf(fid, 'SCALARS vpw_nds float \n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(v(i*df,cont)) '\n']);
    end
      %apw
    fprintf(fid, 'SCALARS apw_nds float \n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
         fprintf(fid,[num2str(a(i*df,cont)) '\n']);
     end
    
    %du
    fprintf(fid, 'VECTORS d_u float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(d(i*df-2,cont)) ' ' num2str(d(i*df-1,cont)) ' 0\n']);
    end
    %Vu
    fprintf(fid, 'VECTORS v_u float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(v(i*df-2,cont)) ' ' num2str(v(i*df-1,cont)) ' 0\n']);
    end
     %Au
     fprintf(fid, 'VECTORS a_u float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
         fprintf(fid,[num2str(a(i*df-2,cont)) ' ' num2str(a(i*df-1,cont)) ' 0\n']);
     end
     %Gamma
     fprintf(fid, 'SCALARS Ep_nds float \n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
         fprintf(fid,[num2str(gamma_nds(i,cont)) '\n']);
     end

    
    fprintf(fid, ['CELL_DATA ' num2str(elements) '\n']);
    %Pressure
    fprintf(fid,'SCALARS Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ps(i*it,cont)) '\n']);
    end
    %Sy
    fprintf(fid,'SCALARS Yield_str float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Sy_tot(i*it,cont)) '\n']);
    end
     %Plastic strain
     fprintf(fid,'SCALARS Ep float\n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:elements
         fprintf(fid,[num2str(gamma(i*it,cont)) '\n']);
     end
    %Pore Pressure
    fprintf(fid,'SCALARS Pore_Pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Pw(i*it,cont)) '\n']);
    end
    
    %Grad Pore Pressure
    fprintf(fid, 'VECTORS Grad_Pore_Pressure float \n');
    for i=1:elements
        fprintf(fid,[num2str(dPw((i*it-1)*sp+1,cont)) ' ' num2str(dPw((i*it-1)*sp+2,cont)) ' 0\n']);
    end
    
    %E11
    fprintf(fid,'SCALARS E11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*it*4-3,cont)) '\n']);
    end
    %E22
    fprintf(fid,'SCALARS E22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*it*4-2,cont)) '\n']);
    end
    %G12
    fprintf(fid,'SCALARS G12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es(i*it*4,cont)) '\n']);
    end
%     %EW11
%     fprintf(fid,'SCALARS EW11 float\n');
%     fprintf(fid, 'LOOKUP_TABLE default\n');
%     for i=1:elements
%         fprintf(fid,[num2str(Es_w(i*5-1,cont)) '\n']);
%     end
%     %EW22
%     fprintf(fid,'SCALARS EW22 float\n');
%     fprintf(fid, 'LOOKUP_TABLE default\n');
%     for i=1:elements
%         fprintf(fid,[num2str(Es_w(i*4,cont)) '\n']);
%     end

    fclose(fid);
end

end
 