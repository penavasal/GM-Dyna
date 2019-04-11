
function read_material

% File: read_material
%   Read the material type and properties from mat.txt
%
% Materials:
% 0-Lin Elastic, 1-Neo-Hookean, 2-VM 3-DP Outer cone, 4-DP Inner cone, 
% 5-DP Plain strain, 6-MCC
%
% Date:
%   Version 1.0   21.03.2018

    global MATERIAL SOLVER GEOMETRY
    
    x_0=GEOMETRY.x_0;
    [~,sp]=size(x_0);
    
    %GEOM
    L=max(x_0(:,1));
    H=max(x_0(:,2));
    
    L0=min(x_0(:,1));
    H0=min(x_0(:,2));
    
    %FILE

    fid = fopen('mat.txt', 'rt'); % opci�n rt para abrir en modo texto
    formato = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s '; % formato de cada l�nea 
    data = textscan(fid, formato, 'HeaderLines', 1);

    a = data{1};
    % Convertir a vector num�rico
    bb = data{2};
    long=length(bb);
    c = data{3};
    
    l=1;
    mats=str2double(bb{l});
    while isnan(mats)
        l=l+1;
        mats=str2double(bb{l});
    end
    MATERIAL.MAT=zeros(31,mats);   % Numero m�ximo de propiedades reconocidas
    MATERIAL.MODEL=zeros(mats,1);
    RANGE=zeros(sp,2*mats); % Range of for materials

    M=0;
    K0=0;
    ks=0;
    kw=0;
    t=l;
    while (M<mats+1) && (t<long)
        t=t+1;
        s1=a{t};
        if strcmp(s1,'MAT')
            M=str2double(bb{t});
            if(M>=mats+1)
                break
            end
            s2=c{t};
            switch s2
                case 'LINEAR_ELASTIC'
                    MATERIAL.MODEL(M)=0;
                    continue
                case 'NEO_HOOKEAN' 
                    MATERIAL.MODEL(M)=1.0;
                    continue
                case 'NEO_HOOKEAN_WRIGGERS' 
                    MATERIAL.MODEL(M)=1.0;
                    continue
                case 'NEO_HOOKEAN_BONET' 
                    MATERIAL.MODEL(M)=1.1;
                    continue
                case 'NEO_HOOKEAN_EHLERS' 
                    MATERIAL.MODEL(M)=1.2;
                    continue
                case 'VON_MISES'
                    MATERIAL.MODEL(M)=2.0;
                    continue
                case 'DRUCKER_PRAGER'
                    MATERIAL.MODEL(M)=2.1;
                    continue
                case 'DRUCKER_PRAGER_O'
                    MATERIAL.MODEL(M)=2.1;
                    continue
                case 'DRUCKER_PRAGER_I'
                    MATERIAL.MODEL(M)=2.2;
                    continue
                case 'DRUCKER_PRAGER_PS'
                    MATERIAL.MATERIAL.MODEL(M)=2.3;
                    continue
                case 'MODIFIED_CAM_CLAY'
                    MATERIAL.MODEL(M)=3.0;
                    continue
                case 'MODIFIED_CAM_CLAY_VISCO'
                    MATERIAL.MODEL(M)=3.1;
                    continue
                otherwise
                    disp('Error, no such material model!')
                    stop
            end
        end
        if strcmp(s1,'X_RANGE')
            d = str2double(bb{t});
            if isnan(d)
                if strcmp(bb{t},'FULL')
                    RANGE(1,M*2-1)=L;
                elseif strcmp(bb{t},'INI')
                    RANGE(1,M*2-1)=L0;
                else
                    disp('Error, wrong material range!')
                    stop
                end
            else
                RANGE(1,M*2-1)=d;
            end
            d = str2double(c{t});
            if isnan(d)
                if strcmp(c{t},'FULL')
                    RANGE(1,M*2)=L;
                elseif strcmp(c{t},'INI')
                    RANGE(1,M*2)=L0;
                else
                    disp('Error, wrong material range!')
                    stop
                end
            else
                RANGE(1,M*2)=d;
            end
            continue
        end
        if strcmp(s1,'Y_RANGE')
            d = str2double(bb{t});
            if isnan(d)
                if strcmp(bb{t},'FULL')
                    RANGE(2,M*2-1)=H;
                elseif strcmp(bb{t},'INI')
                    RANGE(2,M*2-1)=H0;
                else
                    disp('Error, wrong material range!')
                    stop
                end
            else
                RANGE(2,M*2-1)=d;
            end
            d = str2double(c{t});
            if isnan(d)
                if strcmp(c{t},'FULL')
                    RANGE(2,M*2)=H;
                elseif strcmp(c{t},'INI')
                    RANGE(2,M*2)=H0;
                else
                    disp('Error, wrong material range!')
                    stop
                end
            else
                RANGE(2,M*2)=d;
            end
            continue
        end
        if strcmp(s1, '//')
        	continue
        end
        switch s1
            case 'YOUNG'
                MATERIAL.MAT(1,M)=str2double(bb{t});
                continue
            case 'POISSON'
                MATERIAL.MAT(2,M)=str2double(bb{t});
                continue
            case 'DENSITY'
                MATERIAL.MAT(3,M)=str2double(bb{t});
                continue
            case 'SHEAR_MODULUS'
                MATERIAL.MAT(4,M)=str2double(bb{t});
                continue
            case 'LAME_CONSTANT'
                MATERIAL.MAT(5,M)=str2double(bb{t});
                continue
            case 'BULK_MODULUS'
                MATERIAL.MAT(29,M)=str2double(bb{t});
                continue
            case 'CONSTRAINED_MODULUS'
                MATERIAL.MAT(17,M)=str2double(bb{t});
                continue
            case 'WAVE_SPEED'
                MATERIAL.MAT(6,M)=str2double(bb{t});
                continue
            case 'YIELD_STRESS'
                MATERIAL.MAT(7,M)=str2double(bb{t});
                continue
            case 'COHESION'
                MATERIAL.MAT(7,M)=str2double(bb{t});
                continue
            case 'PRECONSOLIDATION'
                MATERIAL.MAT(7,M)=str2double(bb{t});     % Preconsolidation, before OCR
                continue
            case 'INITIAL_PRESSURE'
                MATERIAL.MAT(25,M)=str2double(bb{t});    % Initial pressure
                continue
            case 'HARDENING'
                MATERIAL.MAT(8,M)=str2double(bb{t});
                continue
            case 'HARDENING_EXPONENT'
                MATERIAL.MAT(9,M)=str2double(bb{t});
                continue
            case 'EPSILON0'
                MATERIAL.MAT(10,M)=str2double(bb{t});
                continue
            case 'FRICTION_ANGLE'
                MATERIAL.MAT(11,M)=str2double(bb{t})* pi/180;
                continue
            case 'DILATANCY_ANGLE'
                MATERIAL.MAT(12,M)=str2double(bb{t})* pi/180;
                continue
            case 'VISCOSITY'
                MATERIAL.MAT(13,M)=str2double(bb{t});
                continue
            case 'VISCOSITY_EXPONENT'
                MATERIAL.MAT(14,M)=str2double(bb{t});
                continue
            case 'PERMEABILITY'
                MATERIAL.MAT(15,M)=str2double(bb{t});
                continue
            case 'POROSITY'
                MATERIAL.MAT(16,M)=str2double(bb{t});
                continue
            case 'WATER_BULK_MODULUS'
                MATERIAL.MAT(18,M)=str2double(bb{t});
                continue
            case 'CRITICAL_STATE_LINE'
                MATERIAL.MAT(19,M)=str2double(bb{t});
                continue
            case 'ALPHA_PARAMETER'
                MATERIAL.MAT(20,M)=str2double(bb{t});
                continue
            case 'LAMBDA'
                MATERIAL.MAT(21,M)=str2double(bb{t});
                continue
            case 'KAPPA'
                MATERIAL.MAT(22,M)=str2double(bb{t});
                continue
            case 'INITIAL_VOLUMETRIC_STRAIN'
                MATERIAL.MAT(23,M)=str2double(bb{t});
                continue
            case 'OCR'
                MATERIAL.MAT(24,M)=str2double(bb{t});
                continue
            case 'K0'
                MATERIAL.MAT(26,M)=str2double(bb{t});
                K0=1;
                continue
            case 'KS'
                ks=str2double(bb{t});
                continue
            case 'KW'
                kw=str2double(bb{t});
                continue
            case 'CREEP_INDEX'
                MATERIAL.MAT(30,M)=str2double(bb{t});
                continue
            case 'REFERENCE_TIME'
                MATERIAL.MAT(31,M)=str2double(bb{t});
                continue
            otherwise
                fprintf('Error, no such material property: %s !! \n',s1)
                stop
        end

    end
    
    fclose(fid); 
    
    MATERIAL.number=mats;
    localization(RANGE);
    
    
    for i=1:mats
        if SOLVER.UW
            n=MATERIAL.MAT(16,i);
            if n==0
                fprintf('Error, no porosity !! \n')
            end
            if ks==0
                ks=1e40;       % Pa
                MATERIAL.MAT(27,i)=ks;
            else
                MATERIAL.MAT(27,i)=ks;
            end
            if MATERIAL.MAT(18,i)==0
                if kw==0
                    fprintf('Error, no water bulk modulus !! \n')
                else
                    MATERIAL.MAT(28,i)=kw;
                    MATERIAL.MAT(18,i)=1/(n/kw+(1-n)/ks);
                end
            else
                if kw==0
                    MATERIAL.MAT(28,i)=1/(1/MATERIAL.MAT(18,i)+(1-n)/n/ks);
                else
                    fprintf(...
                    'Error, two different values of the water bulk modulus !! \n')
                end
            end
        end
        if MATERIAL.MODEL(i)<3
            MATERIAL.MAT(:,i)=elastic_tools(MATERIAL.MAT(:,i));
        end
        if MATERIAL.MODEL(i)<3 && MATERIAL.MODEL(i)>=2
            MATERIAL.MAT(:,i)=dp_tools(MATERIAL.MAT(:,i));
        elseif MATERIAL.MODEL(i)<4 && MATERIAL.MODEL(i)>=3
            MATERIAL.MAT(:,i)=mcc_tools(MATERIAL.MAT(:,i),K0);
        end
        
    end




end

function Mat=elastic_tools(Mat)

    global SOLVER
    
    if Mat(4)==0 && Mat(5)==0 && Mat(29)==0
        Mat(4)=Mat(1)/2/(1+Mat(2));             % G 
        Mat(5)=2*Mat(4)*Mat(2)/(1-2*Mat(2));    % Lambda
        Mat(29)=Mat(5)+2*Mat(4);                % K 
    elseif Mat(1)==0 && Mat(5)==0 && Mat(29)==0
        Mat(1)=2*Mat(4)*(1+Mat(2));             % E
        Mat(5)=2*Mat(4)*Mat(2)/(1-2*Mat(2));    % Lambda
        Mat(29)=Mat(5)+2*Mat(4);                % K 
    elseif Mat(1)==0 && Mat(4)==0 && Mat(29)==0
        Mat(1)=Mat(5)*(1+Mat(2))*(1-2*Mat(2))/Mat(2);% E
        Mat(4)=Mat(1)/2/(1+Mat(2));                  % G
        Mat(29)=Mat(5)+2*Mat(4);                     % K
    elseif Mat(1)==0 && Mat(2)==0 && Mat(29)==0
        Mat(1)=Mat(4)*(3*Mat(5)+2*Mat(4))/(Mat(4)+Mat(5)); % E
        Mat(2)=Mat(5)/2/(Mat(5)+Mat(4));                   %nu
        Mat(29)=Mat(5)+2*Mat(4);                           % K
    elseif Mat(1)==0 && Mat(2)==0 && Mat(5)==0
        Mat(5)=Mat(29)-2/3*Mat(4);                         % Lambda
        Mat(1)=Mat(4)*(3*Mat(5)+2*Mat(4))/(Mat(4)+Mat(5)); % E
        Mat(2)=Mat(5)/2/(Mat(5)+Mat(4));                   %nu
    elseif Mat(1)==0 && Mat(4)==0 && Mat(5)==0
        Mat(1)=3*Mat(29)*(1-2*Mat(2));  % E
        Mat(4)=Mat(1)/2/(1+Mat(2));     % G 
        Mat(5)=Mat(29)-2/3*Mat(4);      % Lambda
    elseif Mat(1)==0 && Mat(2)==0 && Mat(4)==0
        Mat(4)=3/2*(Mat(29)-Mat(5));     % G
        Mat(2)=Mat(5)/2/(Mat(5)+Mat(4)); %nu
        Mat(1)=3*Mat(29)*(1-2*Mat(2));   % E
    end
    
    
    % Wave velocity
    Mat(17)=Mat(1)*(1-Mat(2))/((1+Mat(2))*(1-2*Mat(2)));
    if SOLVER.UW
        M=Mat(17)+Mat(18);
    else
        M=Mat(17);
    end
    Mat(6)=sqrt(M/Mat(3));     % cp
    
end

function Mat=dp_tools(Mat)

    if Mat(8)==0
        Mat(8)=1;
        Mat(9)=1e10;
    else
        if Mat(9)==0
            Mat(9)=1;
        end
    end

    Mat(10)= Mat(7)/Mat(8);

end

function Mat=mcc_tools(Mat,K0)

    global SOLVER

    if Mat(19)==0 && Mat(11)==0
        disp('Error, no critical state line!')
        stop
    elseif Mat(19)==0
        Mat(19)=6*sin(Mat(11))/(3-sin(Mat(11)));
    elseif Mat(11)==0
        Mat(11)=asin((3*Mat(19))/(6+Mat(19)));
    end
    
    % OCR
    if Mat(25) && Mat(7)
        Mat(24) = Mat(25)/Mat(7);
    elseif Mat(25) && Mat(24)
        Mat(7) = Mat(25)*Mat(24);
    elseif Mat(7) && Mat(24)
        Mat(25) = Mat(7)/Mat(24);
    elseif Mat(7)
        Mat(25)=Mat(7);
        Mat(24)=1;
    elseif Mat(25)
        Mat(7)=Mat(25);
        Mat(24)=1;
    elseif Mat(24)
        disp('We need more MCC parameters!')
        stop
    end
    
    
    if Mat(21)==0 || Mat(22)==0
        disp('We need more MCC parameters!')
        stop
    end
    
    %%K0
    if Mat(26)==0 && K0==0
        Mat(26)=(1-sin(Mat(11)))*Mat(24)^(sin(Mat(11)));
    end
    
    %%% E
    Mat(1)=2*Mat(4)*(1+Mat(2)); 
    % Wave velocity
    Mat(17)=Mat(1)*(1-Mat(2))/((1+Mat(2))*(1-2*Mat(2)));
    if SOLVER.UW
        M=Mat(17)+Mat(18);
    else
        M=Mat(17);
    end
    Mat(6)=sqrt(M/Mat(3));     % cp
    
end


function localization(R)

    global MATERIAL GEOMETRY
    
    for m=1:MATERIAL.number
        for nodo=1:GEOMETRY.nodes
            tol=GEOMETRY.h_nds(nodo,1)/5;
            if (GEOMETRY.x_0(nodo,2)>=R(2,m*2-1)-tol) ...
                    && (GEOMETRY.x_0(nodo,2)<=R(2,m*2)+tol) 
                if (GEOMETRY.x_0(nodo,1)>=R(1,m*2-1)-tol) ...
                        && (GEOMETRY.x_0(nodo,1)<=R(1,m*2)+tol)
                    MATERIAL.n(nodo)=m;
                end
            end
        end

        for e=1:GEOMETRY.mat_points
            tol=GEOMETRY.h_ini(e,1)/5;
            if (GEOMETRY.xg_0(e,2)>=R(2,m*2-1)-tol) ...
                    && (GEOMETRY.xg_0(e,2)<=R(2,m*2)+tol)
                if (GEOMETRY.xg_0(e,1)>=R(1,m*2-1)-tol) ...
                        && (GEOMETRY.xg_0(e,1)<=R(1,m*2)+tol)
                    MATERIAL.e(e)=m;
                end
            end
        end
        
    end

end



