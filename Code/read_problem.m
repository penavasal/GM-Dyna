
function MAT_POINT=read_problem(str)

% File: read_problem
%   Read some important problem parameters from problem.txt
%
% Date:
%   Version 3.0   26.11.2019

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some possible parameters if they are not read
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global VARIABLE SOLVER TIME GEOMETRY
    
    % GEOMETRY
    filename='data_m';
    filegrid='';
    GRID='';
    pathgeo='Geom';
    DIM=0;
    PLOT_ini=0;
    AMP=1;
        
    SOLVER.REMAPPING = 0;
    SOLVER.LIN = 0;
    SOLVER.AXI=0;
        
    SOLVER.B_BAR=0;
    SOLVER.F_BAR=0;
    SOLVER.F_BAR_W=0;
    SOLVER.Pstab=0;
    
    SOLVER.INITIAL_PORE_PRESSURE=0;
    
    SOLVER.INIT_STEP=0;
    
    SOLVER.TYPE={'' ''};
    
    SOLVER.thickness=1;
    
    SOLVER.FAIL = 0;
    
    SOLVER.BLOCKS = 1;
    
    SOLVER.SMALL = 0;
    
    SOLVER.FRAC = 0;
    
    SOLVER.PHASES = {'' ''};
    
    % VARIABLES
    VARIABLE.g=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read problem.txt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(str, 'rt'); % opción rt para abrir en modo texto
    formato = '%s %s %s %s %s %s %s %s %s'; % formato de cada línea 
    data = textscan(fid, formato, 'HeaderLines', 1);
    a = data{1};
    % Convertir a vector numérico
    b1= cellfun(@str2num, data{2}, 'UniformOutput', false);
    b2= data{2};
    b3= data{3};
    b4= data{4};
    [l,~] = size(b1);
    b=zeros(l,1);
    for i=1:l
        if b1{i}
            b(i)=b1{i};
        else
            b(i)=0;
        end
    end
    t=0;
    while (t<l)
        t=t+1;
        s1=a{t};
        
        if s1=='*'
            if strcmp(b2{t},'TYPE_OF_PROBLEM')
                disp('Reading problem')
                continue
            else
                break
            end
        end
        
        switch s1
            case '//'
                continue
            case 'INIT_FILE'
                val = str2double(b2{t});
                if isnan(val)
                    SOLVER.INIT_file=b2{t};  % Start from a previous stage
                else
                    SOLVER.INIT_file=b(t);
                end
                continue
            case 'INIT_STEP'
                SOLVER.INIT_STEP=b(t);  %
                continue 
            case 'PLOT_INI'
                PLOT_ini=b(t);  % Geometry at the beggining
                continue
            case 'DIMENSION'
                DIM=b(t);        % Flag for the Dimension
                continue
            case 'ELEMENT'
                ELEMENT=b2{t};
                continue
            case 'GRID_TYPE'
                GRID=b2{t};
                continue
            case 'DISCRETIZATION'
                DISC=b(t);  %1- quad (1IP), 2- quad (4IP) 
                            %3- 2P1P0, 4- reverse (3), 5- 4P1P0
                continue
            case 'SCALE'
                AMP=b(t);  % 1, no amplification of the mesh
                continue
            case 'REMAPPING'
                SOLVER.REMAPPING=b(t);
                continue
            case 'FILE'
                filename=b2{t};
                continue
            case 'GRID'
                filegrid=b2{t};
                continue
            case 'PATH_GEOM'
                pathgeo=b2{t};
                continue
            case 'CONFIGURATION'
                if strcmp(b2{t},'PLANE_STRAIN')
                    SOLVER.AXI=0;
                    if SOLVER.thickness==0
                        SOLVER.thickness=1;
                    end
                elseif strcmp(b2{t},'AXISYMMETRIC')
                    SOLVER.AXI=1;
                    SOLVER.thickness=1;
                else
                    disp('Configuration not implemented')
                    stop
                end
                continue
            case 'THICKNESS'
                SOLVER.thickness=b(t);  
                continue
            case 'LINEARIZATION'
                SOLVER.LIN=b(t);
                continue
            case 'FRAMEWORK'
                if strcmp(b2{t},'SMALL_STRAIN')
                    SOLVER.SMALL=1;
                end
                continue
            case 'PROBLEM'
                if strcmp(b2{t},'OTM')
                    SOLVER.TYPE{1}=0;
                elseif strcmp(b2{t},'MPM')
                    SOLVER.TYPE{1}=1;
                elseif strcmp(b2{t},'FE')
                    SOLVER.TYPE{1}=2;
                end
                if b3{t}
                    SOLVER.TYPE{2}=b3{t};
                    if SOLVER.TYPE{2}=='LME'
                        if b4{t}
                           SOLVER.TYPE{3}=b4{t};
                        else
                           fprintf('Check LME file, by Default LME.txt !!\n') 
                           SOLVER.TYPE{3}='LME.txt';
                        end
                    end
                end
                continue
            case 'FORMULATION'
                if strcmp(b2{t},'U')
                    SOLVER.UW=0;
                    SOLVER.PHASES{1,1}='U';
                    SOLVER.PHASES{1,2}=1;
                elseif strcmp(b2{t},'UW')
                    SOLVER.UW=1;
                    SOLVER.PHASES{1,1}='U';
                    SOLVER.PHASES{1,2}=1;
                    SOLVER.PHASES{2,1}='W';
                    SOLVER.PHASES{2,2}=1;
                elseif strcmp(b2{t},'UPw')
                    SOLVER.UW=2;
                    SOLVER.PHASES{1,1}='U';
                    SOLVER.PHASES{1,2}=1;
                    SOLVER.PHASES{2,1}='Pw';
                    SOLVER.PHASES{2,2}=1;
                elseif strcmp(b2{t},'UWPw')
                    SOLVER.UW=3;
                    SOLVER.PHASES{1,1}='U';
                    SOLVER.PHASES{1,2}=1;
                    SOLVER.PHASES{2,1}='W';
                    SOLVER.PHASES{2,2}=1;
                    SOLVER.PHASES{3,1}='Pw';
                    SOLVER.PHASES{3,2}=1;
                end
                if SOLVER.UW>=1
                    if VARIABLE.g==0
                        VARIABLE.g=9.81;         % m/s2
                    end
                end
                continue
            case 'INITIAL_PORE_PRESSURE'
                SOLVER.INITIAL_PORE_PRESSURE=b(t);   
                continue
            case 'GRAVITY'
                VARIABLE.g=b(t);   
                continue
            case 'B_BAR'
                SOLVER.B_BAR=b(t);  
                continue
            case 'F_BAR'
                SOLVER.F_BAR=b(t);  
                continue
            case 'F_BAR_W'
                SOLVER.F_BAR_W=b(t);  
                continue
            case 'PW_STAB'
                SOLVER.Pstab=b(t);  
                continue
            case 'SAVE_FREQUENCY'
                SOLVER.SAVE_I=b(t); % Save info each XX steps
                continue
            case 'FILE_FREQUENCY'
                SOLVER.SAVE_F=b(t); % Save the file each XX steps
                continue
            otherwise
                fprintf('Error, unrecognized parameter: %s !!\n',s1)
                stop
        end
    end
    
    % BLOCKS
    if strcmp(b2{t},'NUMBER_OF_BLOCKS')
        t=t+1;
        SOLVER.BLOCKS=str2double(a{t});
        fprintf('Number of blocks: %i \n',SOLVER.BLOCKS)
    else
        fprintf('Error, unrecognized parameter: %s !!\n',s1)
        stop
    end
    
    % TIME INTEGRATION SCHEME
    TIS=1*ones(SOLVER.BLOCKS,1);
    af=0*ones(SOLVER.BLOCKS,1);
    am=0*ones(SOLVER.BLOCKS,1);
    delta=0*ones(SOLVER.BLOCKS,1);
    alpha=0*ones(SOLVER.BLOCKS,1);
    theta=1*ones(SOLVER.BLOCKS,1);
    rho=0*ones(SOLVER.BLOCKS,1);
    SOLVER.time_factor=1*ones(SOLVER.BLOCKS,1);
    scheme=strings(SOLVER.BLOCKS,1);
    FILES=strings(SOLVER.BLOCKS,4);
    OUT=strings(SOLVER.BLOCKS,1);
    
    SOLVER.OutputType=zeros(1,2);

    while (t<l)
        t=t+1;
        s1=a{t};
        
        if s1=='*'
            if strcmp(b2{t},'BLOCK')
                BLCK=str2double(b3{t});
                fprintf('BLOCK: %i \n',BLCK)
                t=t+1;
                continue
            else
                fprintf('Error, unrecognized parameter: %s !!\n',s1)
                stop
            end
        end
        
        switch s1
            case '//'
                continue
            case 'TIME_FINAL'
                SOLVER.Time_final(BLCK)=b(t);
                continue
            case 'TIME_STEP'
                SOLVER.time_step(BLCK)=b(t);  
                continue
            case 'TIME_FACTOR'
                SOLVER.time_factor(BLCK)=b(t);  
                continue
            case 'DYNAMIC'
                SOLVER.DYN(BLCK)=b(t);  
                if SOLVER.DYN(BLCK)==0
                    scheme(BLCK)='STATIC';
                    TIS(BLCK)=0;
                end
                continue
            case 'SOLVER'
                solver=b2{t};
                if strcmp(solver,'IMPLICIT')
                    SOLVER.IMPLICIT(BLCK)=1;
                elseif strcmp(solver,'EXPLICIT')
                    SOLVER.IMPLICIT(BLCK)=0;
                else
                    fprintf('Error, unrecognized parameter: %s !!\n',s1)
                end
                continue
            case 'OUTPUT'
                OUT(BLCK,1)=b2{t};
                continue
            case 'MATERIAL'
                FILES(BLCK,1)=b2{t};
                continue
            case 'BOUNDARY_CONDITION'
                FILES(BLCK,2)=b2{t};
                continue
            case 'LOAD'
                FILES(BLCK,3)=b2{t};
                continue
            case 'CONTACT'
                FILES(BLCK,4)=b2{t};
                continue
            case 'SCHEME'
                scheme(BLCK)=b2{t};
                switch scheme(BLCK)
                    case 'STATIC'
                        TIS(BLCK)=0;
                        continue
                    case 'NEWMARK1'
                        TIS(BLCK)=1;
                        continue
                    case 'NEWMARK2'
                        TIS(BLCK)=2;
                        continue
                    case 'NEWMARK'
                        TIS(BLCK)=2;
                        continue
                    case 'GENERALIZED_ALPHA'
                        TIS(BLCK)=3;
                        continue
                    case 'HHT'
                        TIS(BLCK)=4;
                        continue
                    case 'WILSON'
                        TIS(BLCK)=5;
                        continue
                    case 'WBZ'
                        TIS(BLCK)=6;
                        continue
                    case 'COLLOCATION'
                        TIS(BLCK)=7;
                        continue
                    case 'NEWMARK_EXPLICIT'
                        TIS(BLCK)=1;
                        continue
                    otherwise
                        disp('Error, no such time integration scheme!')
                        stop
                end
                continue 
            case 'DELTA'
                delta(BLCK)=b(t); 
                continue
            case 'GAMMA'
                delta(BLCK)=b(t); 
                continue
            case 'ALPHA'
                alpha(BLCK)=b(t); 
                continue
            case 'BETA'
                alpha(BLCK)=b(t); 
                continue
            case 'ALPHA_M'
                am(BLCK)=b(t); 
                continue
            case 'ALPHA_F'
                af(BLCK)=b(t); 
                continue
            case 'RHO'
                rho(BLCK)=b(t); 
                continue
            case 'THETA'
                theta(BLCK)=b(t); 
                continue
            case 'NEWTON_RAPHSON_LOOP'
                SOLVER.NR(BLCK)=b(t); %Every NR iterations re-calculate K matrix 
                continue
            case 'NR_TOLERANCE_FORCES'
                SOLVER.r_tolerance(BLCK) =b(t); 
                continue
            case 'NR_TOLERANCE_DISP'
                SOLVER.d_tolerance(BLCK)=b(t);  
                continue
            case 'ITERATIONS' 
                SOLVER.NR_iterations(BLCK)=b(t); 
                continue
            otherwise
                fprintf('Error, unrecognized parameter: %s !!\n',s1)
                stop
        end
    end
    
    fclose(fid); 

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add the geometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [MAT_POINT,NODE_LIST] = read_geometry(...
        ELEMENT,GRID,DIM,PLOT_ini,AMP,filename,filegrid,pathgeo);
    SOLVER.Element=ELEMENT;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add shape function parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [MAT_POINT]=shape_function_calculation(1,MAT_POINT,0,NODE_LIST);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BLOCK problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:SOLVER.BLOCKS
        
        if strcmp(OUT(i),'')
           SOLVER.Output(i)='FILE.mat';
        else
           SOLVER.Output(i)=strcat(OUT(i),'.mat');
        end

        A=exist(SOLVER.Output(i),'file');

        j=0;
        while A==2
            j=j+1;
            SOLVER.Output(i)=strcat(OUT(i),'_',int2str(j),'.mat');
            A=exist(SOLVER.Output(i),'file');
        end
        
        %--------------------------------------------------------------------------
        % PARAMETERS
        %--------------------------------------------------------------------------
        % Read material properties
        read_material(FILES(i,1),i);
         %**************
         % SOLVER
         %**************

         % Solver variables: Time Integration Scheme
         %----------------------------------------------------------------------
         TIME{i}=...
             Time_Scheme(TIS(i),af(i),am(i),delta(i),alpha(i),theta(i),rho(i),i);
 
         % Time-related conditions
         %----------------------------------------------------------------------
         TIME{i}.variables(i);
 
         % BOUNDARY CONDITIONS
         %----------------------------------------------------------
         % FORCES
         %--------------------------------------------------------------------------
         read_load(FILES(i,3),i,NODE_LIST);
         
         % CONSTRAINTS
         %--------------------------------------------------------------------------
         read_boundary(FILES(i,2),i,NODE_LIST);
         
         % Contact conditions  && Potential Surface
         %--------------------------------------------------------------------------
         if isempty(FILES(i,4))
            read_contact; 
         end
    end
end

