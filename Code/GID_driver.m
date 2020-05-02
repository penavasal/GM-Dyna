


function GID_driver(str,str2,steps,rt,folder)

    root=strcat(folder,'/',str,'.mat');
    load(root,'GEOMETRY','SOLVER');
    
    if exist(strcat(folder,'/','GID'),'dir')==0
        mkdir(folder,'GID')
    end

    %[~,NNE]=size(GEOMETRY.elem);
    
    if strcmp(steps,'FULL')
        steps=SOLVER.dim;
    end
    
    print_file_1(str,str2);
    print_file_2(str,str2,steps,rt);
    
    %NNE;



end


function print_file_1(str,str2)

    load(str,'-mat','GEOMETRY');
    [elements,NNE]=size(GEOMETRY.elem);
    elem=GEOMETRY.elem;
    x_a=GEOMETRY.x_0;
    [nodes,sp]=size(x_a);
    
    if NNE==3 || NNE==6
        elemtype='Triangle';
    elseif NNE==4
        elemtype='Quadrilateral';
    elseif NNE==8
        elemtype='Quadrilateral';
    elseif NNE==2
        elemtype='Linear';
    else
        stop;
    end
    
%     type=SOLVER.Element;
%     rr=1;
%     if NNE==4 && strcmp(type,'Q4-4')
%         rr=4;
%     end
    
    flnm=(['GID/' str2 '.post.msh']);
    
    fid = fopen(flnm, 'w'); 
    fprintf(fid, '%25s  \n' , '# GiD Post Mesh for GM-Dyna');
    fprintf(fid,['MESH 1 dimension ' num2str(sp) ' ElemType ' elemtype ...
        ' Nnode ' num2str(NNE) '\n']);
    fprintf(fid,'Coordinates \n');
    for i=1:nodes
        fprintf(fid,[num2str(i,'%10i\t') ' ' num2str(x_a(i,1),'%10.5f\t') ' '...
            num2str(x_a(i,2),'%10.5f\t') ' 0\n']);
    end
    fprintf(fid,'end coordinates \n');
    fprintf(fid,'\n');
    fprintf(fid,'Elements \n');

    for i=1:elements
        fprintf(fid,[num2str(i,'%10i\t') '\t']);
        for j=1:NNE
            fprintf(fid,[num2str(elem(i,j),'%10i\t') '\t']);
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'end elements \n');
end



function print_file_2(str,str2,steps,dd)

    load(str,'-mat','GLOBAL','GEOMETRY','SOLVER');
    [elements,NNE]=size(GEOMETRY.elem);
    
    x_a=GEOMETRY.x_0;
    [nodes,sp]=size(x_a);
    df=GEOMETRY.df;
    
    try
        if strcmp(GEOMETRY.ELEMENT,'Q8P4-4') || strcmp(GEOMETRY.ELEMENT,'Q8P4')
            nodes_p = intersect(GEOMETRY.elem,GEOMETRY.elem_c);
        elseif strcmp(GEOMETRY.ELEMENT,'T6P3-3') || strcmp(GEOMETRY.ELEMENT,'T6P3')
            nodes_p = intersect(GEOMETRY.elem,GEOMETRY.elem_c);
        else
            nodes_p = linspace(1,GEOMETRY.nodes,GEOMETRY.nodes)';
        end
    catch 
       nodes_p = linspace(1,GEOMETRY.nodes,GEOMETRY.nodes)';
    end    
    
    type=SOLVER.Element;
    
    try 
        FRAC=SOLVER.FRAC;
    catch
        FRAC=0;
    end

%     if NNE==3
%         elemtype='Triangle';
%     elseif NNE==4
%         elemtype='Quadrilateral';
%     elseif NNE==2
%         elemtype='Linear';
%     else
%         stop;
%     end
    
    d=GLOBAL.d;
    a=GLOBAL.a;
    v=GLOBAL.v;
    UW=SOLVER.UW;
    
    gamma_nds=GLOBAL.gamma_nds;
    Ps=GLOBAL.Ps;
    Qs=GLOBAL.Qs;
    
    Ss=GLOBAL.Sigma;
    Es=GLOBAL.Es;
    if UW
        Pw=GLOBAL.pw;
        if UW==2
            dPw=GLOBAL.dpw;
        end
    end
    Es_p=GLOBAL.Es_p;
    Sy_tot=GLOBAL.Sy;
    gamma=GLOBAL.gamma;

    if FRAC>0
    	W=GLOBAL.w;
        xi=GLOBAL.status; 
    end

    %Output file name
    flnm=(['GID/' str2 '.post.res']);

    myfile = fopen(flnm,'w');

    fprintf( myfile , '%25s  \n' , 'GiD Post Results File 1.0');
    if NNE==4 || NNE==8
        fprintf( myfile , '%44s  \n' , 'GaussPoints Group1   ElemType  Quadrilateral');
    else
        fprintf( myfile , '%44s  \n' , 'GaussPoints Group1   ElemType  Triangle');
    end
    if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
        fprintf( myfile , '%26s  \n' , 'Number Of Gauss Points:  4');
    elseif strcmp(type,'T3-3') || strcmp(type,'T6P3-3') || strcmp(type,'T6-3') 
        fprintf( myfile , '%26s  \n' , 'Number Of Gauss Points:  3');
    else
        fprintf( myfile , '%26s  \n' , 'Number Of Gauss Points:  1');
    end
    fprintf( myfile , '%29s  \n' , 'Natural Coordinates: Internal');
    fprintf( myfile , '%15s  \n' , 'End GaussPoints');

    if steps==0
        steps=GLOBAL.ste_p-1;
    end

    for cont=1:dd:steps
        % DEGREES of FREEDOM
        % U disp
        fprintf(myfile, ['Result  Displacement SOLID_PHASE ' num2str(GLOBAL.tp(cont)) ' Vector  OnNodes\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:nodes
            fprintf(myfile,[num2str(i,'%10i') ' ' num2str(d((i-1)*df+1,cont)) ' ' num2str(d((i-1)*df+2,cont)) ' 0\n']);
        end
        fprintf( myfile , '%10s  \n' , 'End Values');

        if UW==1
            % W disp
            fprintf(myfile, ['Result  Displacement_water WATER_PHASE ' num2str(GLOBAL.tp(cont)) ' Vector  OnNodes\n']);
            fprintf( myfile , '%6s  \n' ,  'Values');
            for j=1:length(nodes_p)
                i=nodes_p(j);
                fprintf(myfile,[num2str(j,'%10i') ' ' num2str(d((i-1)*df+3,cont)) ' ' num2str(d((i-1)*df+4,cont)) ' 0\n']);
            end
            fprintf( myfile , '%10s  \n' , 'End Values'); 
        end
        if UW==2
            % Pw
            fprintf(myfile, ['Result  Pore_water_Pressure PORE_PRESSURE ' num2str(GLOBAL.tp(cont)) ' Scalar  OnNodes\n']);
            fprintf( myfile , '%6s  \n' ,  'Values');
            for j=1:length(nodes_p)
                i=nodes_p(j);
                fprintf(myfile,[num2str(i,'%10i') ' ' num2str(d((i-1)*df+3,cont)) '\n']);
            end
            fprintf( myfile , '%10s  \n' , 'End Values'); 
        end
        
        % DEGREES of FREEDOM - dor
        % U velo
        fprintf(myfile, ['Result  Velocity SOLID_PHASE ' num2str(GLOBAL.tp(cont)) ' Vector  OnNodes\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:nodes
            fprintf(myfile,[num2str(i,'%10i') ' ' num2str(v((i-1)*df+1,cont)) ' ' num2str(v((i-1)*df+2,cont)) ' 0\n']);
        end
        fprintf( myfile , '%10s  \n' , 'End Values');

        if UW==1
            % W velo
            fprintf(myfile, ['Result  Velocity_water WATER_PHASE ' num2str(GLOBAL.tp(cont)) ' Vector  OnNodes\n']);
            fprintf( myfile , '%6s  \n' ,  'Values');
            for j=1:length(nodes_p)
                i=nodes_p(j);
                fprintf(myfile,[num2str(i,'%10i') ' ' num2str(v((i-1)*df+3,cont)) ' ' num2str(v((i-1)*df+4,cont)) ' 0\n']);
            end
            fprintf( myfile , '%10s  \n' , 'End Values'); 
        end
        if UW==2
            % Pw velo
            fprintf(myfile, ['Result  Pore_water_Pressure_velocity PORE_PRESSURE ' num2str(GLOBAL.tp(cont)) ' Scalar  OnNodes\n']);
            fprintf( myfile , '%6s  \n' ,  'Values');
            for j=1:length(nodes_p)
                i=nodes_p(j);
                fprintf(myfile,[num2str(i,'%10i') ' ' num2str(v((i-1)*df+3,cont)) '\n']);
            end
            fprintf( myfile , '%10s  \n' , 'End Values'); 
        end
        
        % DEGREES of FREEDOM - dor
        % U velo
        fprintf(myfile, ['Result  Acceleration SOLID_PHASE ' num2str(GLOBAL.tp(cont)) ' Vector  OnNodes\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:nodes
            fprintf(myfile,[num2str(i,'%10i') ' ' num2str(a((i-1)*df+1,cont)) ' ' num2str(a((i-1)*df+2,cont)) ' 0\n']);
        end
        fprintf( myfile , '%10s  \n' , 'End Values');

        if UW==1
            % W velo
            fprintf(myfile, ['Result  Acceleration_water WATER_PHASE ' num2str(GLOBAL.tp(cont)) ' Vector  OnNodes\n']);
            fprintf( myfile , '%6s  \n' ,  'Values');
            for j=1:length(nodes_p)
                i=nodes_p(j);
                fprintf(myfile,[num2str(i,'%10i') ' ' num2str(a((i-1)*df+3,cont)) ' ' num2str(a((i-1)*df+4,cont)) ' 0\n']);
            end
            fprintf( myfile , '%10s  \n' , 'End Values'); 
        end
        if UW==2
            % Pw acce
            fprintf(myfile, ['Result  Pore_water_Pressure_acce PORE_PRESSURE ' num2str(GLOBAL.tp(cont)) ' Scalar  OnNodes\n']);
            fprintf( myfile , '%6s  \n' ,  'Values');
            for j=1:length(nodes_p)
                i=nodes_p(j);
                fprintf(myfile,[num2str(i,'%10i') ' ' num2str(a((i-1)*df+3,cont)) '\n']);
            end
            fprintf( myfile , '%10s  \n' , 'End Values'); 
        end   
        % Eq_plastic_strain
        fprintf(myfile, ['Result  Eq_plastic_strain Plastic_Strain ' num2str(GLOBAL.tp(cont)) ' Scalar  OnNodes\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:nodes
            fprintf(myfile,[num2str(i,'%10i') ' ' num2str(gamma_nds(i,cont)) '\n']);
        end
        fprintf( myfile , '%10s  \n' , 'End Values'); 

        % Q
        fprintf(myfile, ['Result  Q Equivalent_stress ' num2str(GLOBAL.tp(cont)) ' Scalar  OnGaussPoints Group1\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:elements
            fprintf(myfile,[num2str(i,'%10i') '\t']);
            if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                for j=1:4
                    if j==1
                        fprintf( myfile ,[num2str(Qs((i-1)*4+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Qs((i-1)*4+j,cont)) '\n']);
                    end
                end
            elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                for j=1:3
                    if j==1
                        fprintf( myfile ,[num2str(Qs((i-1)*3+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Qs((i-1)*3+j,cont)) '\n']);
                    end
                end
            else
                fprintf( myfile ,[num2str(Qs(i,cont)) '\n']);
            end
        end
        fprintf( myfile , '%10s  \n' , 'End Values');
        
        % Ps
        fprintf(myfile, ['Result  P Equivalent_stress ' num2str(GLOBAL.tp(cont)) ' Scalar  OnGaussPoints Group1\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:elements
            fprintf(myfile,[num2str(i,'%10i') '\t']);
            if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                for j=1:4
                    if j==1
                        fprintf( myfile ,[num2str(Ps((i-1)*4+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Ps((i-1)*4+j,cont)) '\n']);
                    end
                end
            elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                for j=1:3
                    if j==1
                        fprintf( myfile ,[num2str(Ps((i-1)*3+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Ps((i-1)*3+j,cont)) '\n']);
                    end
                end
            else
                fprintf( myfile ,[num2str(Ps(i,cont)) '\n']);
            end
        end
        fprintf( myfile , '%10s  \n' , 'End Values');
        
        % Sy
        fprintf(myfile, ['Result  Yield_stress Plasticity ' num2str(GLOBAL.tp(cont)) ' Scalar  OnGaussPoints Group1\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:elements
            fprintf(myfile,[num2str(i,'%10i') '\t']);
            if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                for j=1:4
                    if j==1
                        fprintf( myfile ,[num2str(Sy_tot((i-1)*4+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Sy_tot((i-1)*4+j,cont)) '\n']);
                    end
                end
            elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                for j=1:3
                    if j==1
                        fprintf( myfile ,[num2str(Sy_tot((i-1)*3+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Sy_tot((i-1)*3+j,cont)) '\n']);
                    end
                end
            else
                fprintf( myfile ,[num2str(Sy_tot(i,cont)) '\n']);
            end
        end
        fprintf( myfile , '%10s  \n' , 'End Values');

        % Ep
        fprintf(myfile, ['Result  Equivalent_Plastic_Strain Plasticity ' num2str(GLOBAL.tp(cont)) ' Scalar  OnGaussPoints Group1\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:elements
            fprintf(myfile,[num2str(i,'%10i') '\t']);
            if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                for j=1:4
                    if j==1
                        fprintf( myfile ,[num2str(gamma((i-1)*4+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(gamma((i-1)*4+j,cont)) '\n']);
                    end
                end
            elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                for j=1:3
                    if j==1
                        fprintf( myfile ,[num2str(gamma((i-1)*3+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(gamma((i-1)*3+j,cont)) '\n']);
                    end
                end    
            else
                fprintf( myfile ,[num2str(gamma(i,cont)) '\n']);
            end
        end
        fprintf( myfile , '%10s  \n' , 'End Values');
        
        if FRAC>0
        % Ep
        fprintf(myfile, ['Result  Damage Fracture ' num2str(GLOBAL.tp(cont)) ' Scalar  OnGaussPoints Group1\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:elements
            fprintf(myfile,[num2str(i,'%10i') '\t']);
            if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                for j=1:4
                    if j==1
                        fprintf( myfile ,[num2str(xi((i-1)*4+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(xi((i-1)*4+j,cont)) '\n']);
                    end
                end
            elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                for j=1:3
                    if j==1
                        fprintf( myfile ,[num2str(xi((i-1)*3+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(xi((i-1)*3+j,cont)) '\n']);
                    end
                end    
            else
                fprintf( myfile ,[num2str(xi(i,cont)) '\n']);
            end
        end
        fprintf( myfile , '%10s  \n' , 'End Values');

                % Ep
        fprintf(myfile, ['Result  Strain_energy Fracture ' num2str(GLOBAL.tp(cont)) ' Scalar  OnGaussPoints Group1\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:elements
            fprintf(myfile,[num2str(i,'%10i') '\t']);
            if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                for j=1:4
                    if j==1
                        fprintf( myfile ,[num2str(W((i-1)*4+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(W((i-1)*4+j,cont)) '\n']);
                    end
                end
            elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                for j=1:3
                    if j==1
                        fprintf( myfile ,[num2str(W((i-1)*3+j,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(W((i-1)*3+j,cont)) '\n']);
                    end
                end    
            else
                fprintf( myfile ,[num2str(W(i,cont)) '\n']);
            end
        end
        fprintf( myfile , '%10s  \n' , 'End Values');

        end
        
        if UW
            % Pw
            fprintf(myfile, ['Result  Pore_Pressure PORE_PRESSURE ' num2str(GLOBAL.tp(cont)) ' Scalar  OnGaussPoints Group1\n']);
            fprintf( myfile , '%6s  \n' ,  'Values');
            for i=1:elements
                fprintf(myfile,[num2str(i,'%10i') '\t']);
                if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                    for j=1:4
                        if j==1
                            fprintf( myfile ,[num2str(Pw((i-1)*4+j,cont)) '\n']);
                        else
                            fprintf( myfile ,['\t' num2str(Pw((i-1)*4+j,cont)) '\n']);
                        end
                    end
                elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                    for j=1:3
                        if j==1
                            fprintf( myfile ,[num2str(Pw((i-1)*3+j,cont)) '\n']);
                        else
                            fprintf( myfile ,['\t' num2str(Pw((i-1)*3+j,cont)) '\n']);
                        end
                    end 
                else
                    fprintf( myfile ,[num2str(Pw(i,cont)) '\n']);
                end
            end
            fprintf( myfile , '%10s  \n' , 'End Values');
            if UW==2
                % Pw
                fprintf(myfile, ['Result  Pore_Pressure_Gradient PORE_PRESSURE ' num2str(GLOBAL.tp(cont)) ' Vector  OnGaussPoints Group1\n']);
                fprintf( myfile , '%6s  \n' ,  'Values');
                for i=1:elements
                    fprintf(myfile,[num2str(i,'%10i') '\t']);
                    if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                        for j=1:4
                            if j==1
                                fprintf( myfile ,[num2str(dPw(((i-1)*4+j-1)*sp+1,cont)) ' ' num2str(dPw(((i-1)*4+j-1)*sp+2,cont)) ' 0\n']);
                            else
                                fprintf( myfile ,['\t' num2str(dPw(((i-1)*4+j-1)*sp+1,cont)) ' ' num2str(dPw(((i-1)*4+j-1)*sp+2,cont)) ' 0\n']);
                            end
                        end
                    elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                        for j=1:3
                            if j==1
                                fprintf( myfile ,[num2str(dPw(((i-1)*3+j-1)*sp+1,cont)) ' ' num2str(dPw(((i-1)*3+j-1)*sp+2,cont)) ' 0\n']);
                            else
                                fprintf( myfile ,['\t' num2str(dPw(((i-1)*3+j-1)*sp+1,cont)) ' ' num2str(dPw(((i-1)*3+j-1)*sp+2,cont)) ' 0\n']);
                            end
                        end 
                    else
                        fprintf( myfile ,[num2str(dPw((i-1)*sp+1,cont)) ' ' num2str(dPw((i-1)*sp+2,cont)) ' 0\n']);
                    end
                end
                fprintf( myfile , '%10s  \n' , 'End Values');
            end
        end
        
        % Stress
        fprintf(myfile, ['Result  Effective_Stress Stress ' num2str(GLOBAL.tp(cont)) ' PlainDeformationMatrix  OnGaussPoints Group1\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:elements
            fprintf(myfile,[num2str(i,'%10i') '\t']);
            if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                for j=1:4
                    if j==1
                        fprintf( myfile ,[num2str(Ss((((i-1)*4+j)-1)*4+1,cont)) ' ' num2str(Ss((((i-1)*4+j)-1)*4+2,cont))...
                            ' ' num2str(Ss((((i-1)*4+j)-1)*4+4,cont)) ' ' num2str(Ss((((i-1)*4+j)-1)*4+3,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Ss((((i-1)*4+j)-1)*4+1,cont)) ' ' num2str(Ss((((i-1)*4+j)-1)*4+2,cont))...
                            ' ' num2str(Ss((((i-1)*4+j)-1)*4+4,cont)) ' ' num2str(Ss((((i-1)*4+j)-1)*4+3,cont)) '\n']);
                    end
                end
            elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                for j=1:3
                    if j==1
                        fprintf( myfile ,[num2str(Ss((((i-1)*3+j)-1)*4+1,cont)) ' ' num2str(Ss((((i-1)*3+j)-1)*4+2,cont))...
                            ' ' num2str(Ss((((i-1)*3+j)-1)*4+4,cont)) ' ' num2str(Ss((((i-1)*3+j)-1)*4+3,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Ss((((i-1)*3+j)-1)*4+1,cont)) ' ' num2str(Ss((((i-1)*3+j)-1)*4+2,cont))...
                            ' ' num2str(Ss((((i-1)*3+j)-1)*4+4,cont)) ' ' num2str(Ss((((i-1)*3+j)-1)*4+3,cont)) '\n']);
                    end
                end 
            else
                fprintf( myfile ,[num2str(Ss((i-1)*4+1,cont)) ' ' num2str(Ss((i-1)*4+2,cont))...
                    ' ' num2str(Ss((i-1)*4+3,cont)) ' ' num2str(Ss((i-1)*4+4,cont)) '\n']);
            end
        end
        fprintf( myfile , '%10s  \n' , 'End Values');
        
        % Strain
        fprintf(myfile, ['Result  Total_Strain Strain ' num2str(GLOBAL.tp(cont)) ' PlainDeformationMatrix  OnGaussPoints Group1\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:elements
            fprintf(myfile,[num2str(i,'%10i') '\t']);
            if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                for j=1:4
                    if j==1
                        fprintf( myfile ,[num2str(Es((((i-1)*4+j)-1)*4+1,cont)) ' ' num2str(Es((((i-1)*4+j)-1)*4+2,cont))...
                            ' ' num2str(Es((((i-1)*4+j)-1)*4+4,cont)) ' ' num2str(Es((((i-1)*4+j)-1)*4+3,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Es((((i-1)*4+j)-1)*4+1,cont)) ' ' num2str(Es((((i-1)*4+j)-1)*4+2,cont))...
                            ' ' num2str(Es((((i-1)*4+j)-1)*4+4,cont)) ' ' num2str(Es((((i-1)*4+j)-1)*4+3,cont)) '\n']);
                    end
                end
            elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                for j=1:3
                    if j==1
                        fprintf( myfile ,[num2str(Es((((i-1)*3+j)-1)*4+1,cont)) ' ' num2str(Es((((i-1)*3+j)-1)*4+2,cont))...
                            ' ' num2str(Es((((i-1)*3+j)-1)*4+4,cont)) ' ' num2str(Es((((i-1)*3+j)-1)*4+3,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Es((((i-1)*3+j)-1)*4+1,cont)) ' ' num2str(Es((((i-1)*3+j)-1)*4+2,cont))...
                            ' ' num2str(Es((((i-1)*3+j)-1)*4+4,cont)) ' ' num2str(Es((((i-1)*3+j)-1)*4+3,cont)) '\n']);
                    end
                end
            else
                fprintf( myfile ,[num2str(Es((i-1)*4+1,cont)) ' ' num2str(Es((i-1)*4+2,cont))...
                    ' ' num2str(Es((i-1)*4+4,cont)) ' ' num2str(Es((i-1)*4+3,cont)) '\n']);
            end
        end
        fprintf( myfile , '%10s  \n' , 'End Values');
        
        % Plastic_Strain
        fprintf(myfile, ['Result  Plastic_Strain Strain ' num2str(GLOBAL.tp(cont)) ' PlainDeformationMatrix  OnGaussPoints Group1\n']);
        fprintf( myfile , '%6s  \n' ,  'Values');
        for i=1:elements
            fprintf(myfile,[num2str(i,'%10i') '\t']);
            if strcmp(type,'Q4-4') || strcmp(type,'Q8-4') || strcmp(type,'Q8P4-4') 
                for j=1:4
                    if j==1
                        fprintf( myfile ,[num2str(Es_p((((i-1)*4+j)-1)*4+1,cont)) ' ' num2str(Es_p((((i-1)*4+j)-1)*4+2,cont))...
                            ' ' num2str(Es_p((((i-1)*4+j)-1)*4+4,cont)) ' ' num2str(Es_p((((i-1)*4+j)-1)*4+3,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Es_p((((i-1)*4+j)-1)*4+1,cont)) ' ' num2str(Es_p((((i-1)*4+j)-1)*4+2,cont))...
                            ' ' num2str(Es_p((((i-1)*4+j)-1)*4+4,cont)) ' ' num2str(Es_p((((i-1)*4+j)-1)*4+3,cont)) '\n']);
                    end
                end
            elseif strcmp(type,'T3-3') || strcmp(type,'T6-3') || strcmp(type,'T6P3-3') 
                for j=1:3
                    if j==1
                        fprintf( myfile ,[num2str(Es_p((((i-1)*3+j)-1)*4+1,cont)) ' ' num2str(Es_p((((i-1)*3+j)-1)*4+2,cont))...
                            ' ' num2str(Es_p((((i-1)*3+j)-1)*4+4,cont)) ' ' num2str(Es_p((((i-1)*3+j)-1)*4+3,cont)) '\n']);
                    else
                        fprintf( myfile ,['\t' num2str(Es_p((((i-1)*3+j)-1)*4+1,cont)) ' ' num2str(Es_p((((i-1)*3+j)-1)*4+2,cont))...
                            ' ' num2str(Es_p((((i-1)*3+j)-1)*4+4,cont)) ' ' num2str(Es_p((((i-1)*3+j)-1)*4+3,cont)) '\n']);
                    end
                end
            else
                fprintf( myfile ,[num2str(Es_p((i-1)*4+1,cont)) ' ' num2str(Es_p((i-1)*4+2,cont))...
                    ' ' num2str(Es_p((i-1)*4+3,cont)) ' ' num2str(Es_p((i-1)*4+4,cont)) '\n']);
            end
        end
        fprintf( myfile , '%10s  \n' , 'End Values');
        
    end
    fclose( myfile );
end

