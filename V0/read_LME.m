function read_LME

% File: read_LME
%   Read some important parameters of the shape function from LME.txt
%
% Date:
%   Version 1.0   16.05.2018

    global LME_param


    % As described in [1], gamma measures the degree of locality, the larger,
    % the closer to Delaunay finite elements

    gamma_lme=3.5;
    gamma_top=0.9;
        
    % This sets the numerical threshold for the support of the shape functions
    target_zero=1.e-6; 
    % This is the tolerance for the Newton iteration to compute the shape
    % functions
    TolLag=max(2*eps,1.e-18);
    % Nelder Mead method to obtain minimum lambda in shape functions
    Nelder=1;
    
    %Tolerance to make new re-mapping or not and proportion to do it
    tol_search = 0.5;       % Optimum [0.4-0.7]
    prop = 0.1;             % Proportion of gamma decreasing
    
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Read LME.txt
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen('LME.txt', 'rt'); % opción rt para abrir en modo texto
    formato = '%s %s %s %s %s %s %s %s %s'; % formato de cada línea 
    data = textscan(fid, formato, 'HeaderLines', 1);
    a = data{1};
    % Convertir a vector numérico
    b1= cellfun(@str2num, data{2}, 'UniformOutput', false);
    b2= data{2};
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
            case 'GAMMA_LME'
                gamma_lme=b(t); 
                continue
            case 'GAMMA_TOP'
                gamma_top=b(t);
                continue
            case 'TARGET_ZERO'
                target_zero=b(t);
                continue
            case 'TOL_LAG'
                TolLag=max(b(t),TolLag);  
                continue
            case 'TOL_SEARCH'
                tol_search=b(t);  
                continue
            case 'WRAPPER'
                if strcmp(b2{t},'NELDER') || strcmp(b2{t},'NELDER_MEAD')
                    Nelder=1;
                elseif strcmp(b2{t},'NR') || strcmp(b2{t},'NEWTON_RAPHSON')
                    Nelder=0;
                end
                continue
            case 'PROPORTION'
                prop=b(t);  
                continue
            otherwise
                fprintf('Error, unrecognized parameter: %s !!\n',s1)
                stop
        end
    end
    
    fclose(fid); 
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Save LME_param
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    LME_param(1)=gamma_lme;
    LME_param(2)=target_zero;
    LME_param(3)=max(2*eps,TolLag);
    LME_param(4)=Nelder;
    LME_param(5)=tol_search;
    LME_param(6)=prop;
    LME_param(7)=gamma_top;

end