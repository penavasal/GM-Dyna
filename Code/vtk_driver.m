function vtk_driver(str,str2,steps,h,rt,folder)

    root=strcat(folder,'/',str,'.mat');

    load(root,'SOLVER');
   
    if strcmp(steps,'FULL')
        steps=SOLVER.dim;
    end
    
    vtk_writter(str,str2,steps,h,rt);

end

function vtk_writter(str,str2,steps,sc,dd)

%clear
load(str,'-mat','GLOBAL','GEOMETRY','SOLVER');
if strcmp(GEOMETRY.ELEMENT,'Q8P4-4') || strcmp(GEOMETRY.ELEMENT,'Q4-4') ...
        || strcmp(GEOMETRY.ELEMENT,'Q8-4')
    it=4;
elseif strcmp(GEOMETRY.ELEMENT,'T6P33') || strcmp(GEOMETRY.ELEMENT,'T3-3') ...
        || strcmp(GEOMETRY.ELEMENT,'T6-3')
        it=3;
else
    it=1;
end
[elements,NNE]=size(GEOMETRY.elem);
x_a=GEOMETRY.x_0;
sp=GEOMETRY.sp;
df=GEOMETRY.df;
%[ste_p,~]=size(tp);

    if NNE==4 || NNE==3
        elem=GEOMETRY.elem;
        elem2=GEOMETRY.elem;
        nodes_p = linspace(1,GEOMETRY.nodes,GEOMETRY.nodes)';
        if NNE==4
            elid=' 9 \n';
        else
            elid=' 5 \n';
        end
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
        elid=' 9 \n';
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
        elid=' 5 \n';
    end
    Ss=GLOBAL.Sigma;
    Es=GLOBAL.Es;
    Ps=GLOBAL.Ps;
    Qs=GLOBAL.Qs;
    Es_p=GLOBAL.Es_p;
    Sy_tot=GLOBAL.Sy;
    gamma_nds=GLOBAL.gamma_nds;
    gamma=GLOBAL.gamma;
    d=GLOBAL.d;
    a=GLOBAL.a;
    v=GLOBAL.v;
    
    if SOLVER.UW==1
        Pw=GLOBAL.pw;
        if SOLVER.SMALL==1
            Es_w=GLOBAL.Es_w;
        end
    elseif SOLVER.UW==2 || SOLVER.UW==3
        Pw=GLOBAL.pw;
        dPw=GLOBAL.dpw;
    end
    
    if SOLVER.FRAC
        status=GLOBAL.status;
        W=GLOBAL.w;
    end
    
    if steps==0
        steps=GLOBAL.ste_p-1;
    end

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
        fprintf(fid,[num2str(x_a(i,1)+sc*d((i-1)*df+1,cont)) ' ' num2str(x_a(i,2)+sc*d((i-1)*df+2,cont)) ' 0\n']);
    end
    %grid	
    fprintf(fid, ['CELLS ' num2str(elements) ' ' num2str(elements*(NNE+1)) '\n']);
    for i=1:elements
        if NNE==4
            fprintf(fid,[num2str(NNE) ' ' num2str(elem2(i,1)-1) ' ' num2str(elem2(i,2)-1) ' ' num2str(elem2(i,3)-1) ' ' num2str(elem2(i,4)-1) '\n']);
        elseif NNE==3
            fprintf(fid,[num2str(NNE) ' ' num2str(elem(i,1)-1) ' ' num2str(elem(i,2)-1) ' ' num2str(elem(i,3)-1) '\n']);
        end
    end
    fprintf(fid, ['CELL_TYPES ' num2str(elements) '\n']);
    for i=1:elements
        fprintf(fid,elid);
    end
    fprintf(fid, ['POINT_DATA ' num2str(length(nodes_p)) '\n']);
    
    %du
    fprintf(fid, 'VECTORS d_u float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(d((i-1)*df+1,cont)) ' ' num2str(d((i-1)*df+2,cont)) ' 0\n']);
    end
    %Vu
    fprintf(fid, 'VECTORS v_u float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
        fprintf(fid,[num2str(v((i-1)*df+1,cont)) ' ' num2str(v((i-1)*df+2,cont)) ' 0\n']);
    end
    %Au
    fprintf(fid, 'VECTORS a_u float \n');
    for j=1:length(nodes_p)
        i=nodes_p(j);
         fprintf(fid,[num2str(a((i-1)*df+1,cont)) ' ' num2str(a((i-1)*df+2,cont)) ' 0\n']);
    end
    
    if SOLVER.UW==1 || SOLVER.UW==3
        %dw
        fprintf(fid, 'VECTORS d_w float \n');
        for j=1:length(nodes_p)
            i=nodes_p(j);
            fprintf(fid,[num2str(d((i-1)*df+1+sp,cont)) ' ' num2str(d((i-1)*df+2+sp,cont)) ' 0\n']);
        end
        %Vw
        fprintf(fid, 'VECTORS v_w float \n');
        for j=1:length(nodes_p)
            i=nodes_p(j);
            fprintf(fid,[num2str(v((i-1)*df+1+sp,cont)) ' ' num2str(v((i-1)*df+2+sp,cont)) ' 0\n']);
        end
        %Aw
         fprintf(fid, 'VECTORS a_w float \n');
        for j=1:length(nodes_p)
            i=nodes_p(j);
             fprintf(fid,[num2str(a((i-1)*df+1+sp,cont)) ' ' num2str(a((i-1)*df+2+sp,cont)) ' 0\n']);
        end
    elseif SOLVER.UW==2 || SOLVER.UW==3
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
    %Q
    fprintf(fid,'SCALARS Q float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Qs(i*it,cont)) '\n']);
    end
    %Sy
    fprintf(fid,'SCALARS Yield_str float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Sy_tot(i*it,cont)) '\n']);
    end
    
    if SOLVER.FRAC>0
        %Status
        fprintf(fid,'SCALARS Status float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(status(i*it,cont)) '\n']);
        end
        %W
        fprintf(fid,'SCALARS W float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(W(i*it,cont)) '\n']);
        end
    end
    
     %Plastic strain
     fprintf(fid,'SCALARS Ep float\n');
     fprintf(fid, 'LOOKUP_TABLE default\n');
     for i=1:elements
         fprintf(fid,[num2str(gamma(i*it,cont)) '\n']);
     end
     if SOLVER.UW>0
        %Pore Pressure
        fprintf(fid,'SCALARS Pore_Pressure float\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i=1:elements
            fprintf(fid,[num2str(Pw(i*it,cont)) '\n']);
        end
        if SOLVER.UW==2 || SOLVER.UW==3
                %Grad Pore Pressure
            fprintf(fid, 'VECTORS Grad_Pore_Pressure float \n');
            for i=1:elements
                fprintf(fid,[num2str(dPw((i*it-1)*sp+1,cont)) ' ' num2str(dPw((i*it-1)*sp+2,cont)) ' 0\n']);
            end
        elseif SOLVER.UW==1 && SOLVER.SMALL==1
            %EW11
            fprintf(fid,'SCALARS EW11 float\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            for i=1:elements
                fprintf(fid,[num2str(Es_w(i*it*4-3,cont)) '\n']);
            end
            %EW22
            fprintf(fid,'SCALARS EW22 float\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            for i=1:elements
                fprintf(fid,[num2str(Es_w(i*it*4-2,cont)) '\n']);
            end
        end
        
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
    
        %E11_p
    fprintf(fid,'SCALARS E11_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*it*4-3,cont)) '\n']);
    end
    %E22_p
    fprintf(fid,'SCALARS E22_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*it*4-2,cont)) '\n']);
    end
    %E33_p
    fprintf(fid,'SCALARS E33_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*it*4-1,cont)) '\n']);
    end
    %G12_p
    fprintf(fid,'SCALARS G12_p float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Es_p(i*it*4,cont)) '\n']);
    end
    
    %S11
    fprintf(fid,'SCALARS S11 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*it*4-3,cont)) '\n']);
    end
    %S22
    fprintf(fid,'SCALARS S22 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*it*4-2,cont)) '\n']);
    end
    %S33
    fprintf(fid,'SCALARS S33 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*it*4-1,cont)) '\n']);
    end
    %T12
    fprintf(fid,'SCALARS T12 float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i=1:elements
        fprintf(fid,[num2str(Ss(i*it*4,cont)) '\n']);
    end



    fclose(fid);
end

end
