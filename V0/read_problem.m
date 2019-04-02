
function read_problem

% File: read_problem
%   Read some important problem parameters from problem.txt
%
% Date:
%   Version 2.0   02.04.2019

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some possible parameters if they are not read
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global VARIABLE SOLVER TI_param
    
    % GEOMETRY
    filename='data_m';
    pathgeo='Geom';
    DIM=0;
    PLOT_ini=0;
    AMP=1;
    
    % TIME INTEGRATION SCHEME
    TIS=1;
    af=0;
    am=0;
    delta=0;
    alpha=0;
    theta=1;
    
    SOLVER.REMAPPING = 0;
    SOLVER.LIN = 0;
    SOLVER.AXI=0;
    
    ELEMENT={'',''};
    
    SOLVER.time_factor=1;
    
    SOLVER.B_BAR=0;
    SOLVER.F_BAR=0;
    SOLVER.F_BAR_W=0;
    
    SOLVER.INITIAL_COND=zeros(2,1);
    
    SOLVER.INIT_STEP=0;
    
    SOLVER.thickness=1;
    
    % VARIABLES
    VARIABLE.rho_w=0;
    VARIABLE.g=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read problem.txt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen('problem.txt', 'rt'); % opción rt para abrir en modo texto
    formato = '%s %s %s %s %s %s %s %s %s'; % formato de cada línea 
    data = textscan(fid, formato, 'HeaderLines', 1);
    a = data{1};
    % Convertir a vector numérico
    b1= cellfun(@str2num, data{2}, 'UniformOutput', false);
    b2= data{2};
    b3= data{3};
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
                ELEMENT{1}=b2{t};
                if b3{t}
                    ELEMENT{2}=b3{t};    
                end
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
                f2=1;
                continue
            case 'PATH_GEOM'
                pathgeo=b2{t};
                f1=1;
                continue
            case 'CONFIGURATION'
                if strcmp(b2{t},'PLAIN_STRAIN')
                    if SOLVER.thickness==0
                        SOLVER.thickness=1;
                    end
                elseif strcmp(b2{t},'AXISYMMETRIC')
                    SOLVER.AXI=1;
                    SOLVER.thickness=1;
                end
                continue
            case 'THICKNESS'
                SOLVER.thickness=b(t);  
                continue
            case 'LINEARIZATION'
                SOLVER.LIN=b(t);   
                continue
            case 'PROBLEM'
                if strcmp(b2{t},'OTM')
                    SOLVER.TYPE=0;
                elseif strcmp(b2{t},'MPM')
                    SOLVER.TYPE=1;
                elseif strcmp(b2{t},'FEM')
                    SOLVER.TYPE=2;
                    fprintf('Error, NOT IMPLEMENTED FEM YET!!\n')
                    stop
                end
                continue
            case 'FORMULATION'
                if strcmp(b2{t},'U')
                    SOLVER.UW=0;
                elseif strcmp(b2{t},'UW')
                    SOLVER.UW=1;
                elseif strcmp(b2{t},'UPw')
                    SOLVER.UW=2;
                elseif strcmp(b2{t},'UWPw')
                    SOLVER.UW=3;
                end
                if SOLVER.UW>=1
                    if VARIABLE.rho_w==0
                        VARIABLE.rho_w=1000;
                    end
                    if VARIABLE.g==0
                        VARIABLE.g=9.81;         % m/s2
                    end
                end
                continue
            case 'WATER_DENSITY'
                VARIABLE.rho_w=b(t);
                continue
            case 'INITIAL_GRAVITY'
                if b(t)==1
                    SOLVER.INITIAL_COND(1)=1;   %Multiplies g
                    if VARIABLE.g==0
                        VARIABLE.g=9.81;         % m/s2
                    end
                elseif b(t)==0
                    SOLVER.INITIAL_COND(1)=0;   %Not g at the beginning
                    if VARIABLE.g==0
                        VARIABLE.g=9.81;         % m/s2
                    end
                end     
                continue
            case 'INITIAL_PORE_PRESSURE'
                SOLVER.INITIAL_COND(2)=b(t);   
                continue
            case 'INITIAL_DISPLACEMENT'
                if strcmp(b2{t},'YES')
                    SOLVER.INITIAL_d=1;
                else
                    SOLVER.INITIAL_d=0;
                end
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
            case 'TIME_FINAL'
                SOLVER.Time_final=b(t);
                continue
            case 'TIME_STEP'
                SOLVER.time_step=b(t);  
                continue
            case 'TIME_FACTOR'
                SOLVER.time_factor=b(t);  
                continue
            case 'SOLVER'
                solver=b2{t};
                if strcmp(solver,'IMPLICIT')
                    SOLVER.IMPLICIT=1;
                elseif strcmp(solver,'EXPLICIT')
                    SOLVER.IMPLICIT=0;
                else
                    fprintf('Error, unrecognized parameter: %s !!\n',s1)
                end
                continue
            case 'SAVE_FREQUENCY'
                SOLVER.SAVE_I=b(t); % Save info each XX steps
                continue
            case 'FILE_FREQUENCY'
                SOLVER.SAVE_F=b(t); % Save the file each XX steps
                continue
            case 'SCHEME'
                scheme=b2{t};
                switch scheme
                    case 'NEWMARK1'
                        TIS=1;
                        continue
                    case 'NEWMARK2'
                        TIS=2;
                        continue
                    case 'GENERALIZED_ALPHA'
                        TIS=3;
                        continue
                    case 'HHT'
                        TIS=4;
                        continue
                    case 'WILSON'
                        TIS=5;
                        continue
                    case 'WBZ'
                        TIS=6;
                        continue
                    case 'COLLOCATION'
                        TIS=7;
                        continue
                    case 'NEWMARK_EXPLICIT'
                        delta=0.5;     %gamma
                        TIS=1;
                        continue
                    otherwise
                        disp('Error, no such time integration scheme!')
                        stop
                end
                continue 
            case 'DELTA'
                delta=b(t); 
                continue
            case 'ALPHA'
                alpha=b(t); 
                continue
            case 'ALPHA_M'
                am=b(t); 
                continue
            case 'ALPHA_F'
                af=b(t); 
                continue
            case 'RHO'
                rho=b(t); 
                continue
            case 'THETA'
                theta=b(t); 
                continue
            case 'NEWTON_RAPHSON_LOOP'
                SOLVER.NR=b(t); %Every NR iterations re-calculate K matrix 
                continue
            case 'NR_TOLERANCE_RELATIVE'
                SOLVER.rel_tolerance =b(t); 
                continue
            case 'NR_TOLERANCE_ABSOLUTE'
                SOLVER.abs_tolerance=b(t);  
                continue
            case 'ITERATIONS' 
                SOLVER.NR_iterations=b(t); 
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
%     if f2==0
%         filename='data_m';
%     end
    
    read_geometry(SOLVER.UW,ELEMENT,DIM,PLOT_ini,AMP,filename,pathgeo);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%SOLVER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SOLVER.IMPLICIT==0 && TIS~=1
        disp('error with the explicit time integration scheme,')
        disp('changed to Newmark');
        TIS=1;
    end
    %----------------------------------------------------------------------
    % Solver variables: Time Integration Scheme
    %----------------------------------------------------------------------
    if TIS==1
        alpha=0;
    elseif TIS==3 || TIS==4 || TIS==6
        if rho
            if am || af
                disp('Error, do I take alpha_f, alpha_m or rho??')
            else
                 if TIS==3                   %GENERALIZED ALPHA
                    af=rho/(1+rho);
                    am=(2*rho-1)/(1+rho);

                    delta=0.5+af-am;
                    alpha=0.25*(1-am+af)^2;
                elseif TIS==4                  %HHT
                    am=0;
                    af=(1-rho)/(1+rho);

                    delta=(1+2*af)/2;
                    alpha=0.25*(1+af)^2;      

                elseif TIS==6                   %WBZ
                    af=0;
                    am=(rho-1)/(1+rho);

                    delta=0.5-am;
                    alpha=0.25*(1-am)^2;
                 end
            end
        elseif am || af
             if TIS==3                   %GENERALIZED ALPHA
                delta=0.5+af-am;
                alpha=0.25*(1-am+af)^2;
            elseif TIS==4                  %HHT
                delta=(1+2*af)/2;
                alpha=0.25*(1+af)^2;      
            elseif TIS==6                   %WBZ
                delta=0.5-am;
                alpha=0.25*(1-am)^2;
             end
        end
    elseif TIS==5                  %WILSON-THETA
        if theta==0
            disp('Error, theta cannot be zero for Wilson')
        end
        alpha=1/4;
        delta=1/2;
    elseif TIS==7                   %COLLOCATION METHOD
        if theta==0 
            disp('Error, theta cannot be zero for Collocation')
        end
        if alpha==0 
            disp('Error, alpha cannot be zero for Collocation')
        end
        if delta==0 
            disp('Error, delta cannot be zero for Collocation')
        end
    end

    TI_param=[af,am,delta,alpha,theta];
    
end
