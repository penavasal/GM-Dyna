
function read_material(filetxt,BLCK)

% File: read_material
%   Read the material type and properties from mat.txt
%
% Materials:
% 0-Lin Elastic, 1-Neo-Hookean, 2-VM 3-DP Outer cone, 4-DP Inner cone, 
% 5-DP Plain strain, 6-MCC
%
% Date:
%   Version 2.0   25.11.2019

    global MATERIAL SOLVER GEOMETRY
    
    x_0=GEOMETRY.x_0;
    [~,sp]=size(x_0);
    
    mat_e=GEOMETRY.material;
    GEOMETRY.body=zeros(length(mat_e),1);
    
    %GEOM
    L=max(x_0(:,1));
    H=max(x_0(:,2));
    
    L0=min(x_0(:,1));
    H0=min(x_0(:,2));
    
    %FILE

    fid = fopen(filetxt, 'rt'); % opci�n rt para abrir en modo texto
    formato = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s '; % formato de cada l�nea 
    data = textscan(fid, formato, 'HeaderLines', 1);

    a = data{1};
    % Convertir a vector numérico
    bb = data{2};
    bb2=str2double(bb);
    long=length(bb);
    c = data{3};
    
    l=1;
    mats=str2double(bb{l});
    while isnan(mats)
        l=l+1;
        mats=str2double(bb{l});
    end
    
    mats2=max(mat_e);
    if mats2<=mats
        MATERIAL(BLCK).number=mats;
    else
        fprintf('Error, need more materials \n')
        stop
    end
    
    BODIES(mats)=0;
    
    MATERIAL(BLCK).MAT=cell(64,mats);   % Numero m¿ximo de propiedades reconocidas
    MATERIAL(BLCK).MODEL=zeros(mats,2);
    %RANGE=zeros(sp,2*mats); % Range of for materials
    
    MATERIAL(BLCK).MAT(61,:)={1};

    M=0;
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
                    MATERIAL(BLCK).MODEL(M)=0;
                    continue
                case 'NEO_HOOKEAN' 
                    MATERIAL(BLCK).MODEL(M)=1.0;
                    continue
                case 'NEO_HOOKEAN_WRIGGERS' 
                    MATERIAL(BLCK).MODEL(M)=1.0;
                    continue
                case 'NEO_HOOKEAN_BONET' 
                    MATERIAL(BLCK).MODEL(M)=1.1;
                    continue
                case 'NEO_HOOKEAN_EHLERS' 
                    MATERIAL(BLCK).MODEL(M)=1.2;
                    continue
                case 'VON_MISES'
                    MATERIAL(BLCK).MODEL(M)=2.0;
                    continue
                case 'DRUCKER_PRAGER'
                    MATERIAL(BLCK).MODEL(M)=2.1;
                    continue
                case 'DRUCKER_PRAGER_O'
                    MATERIAL(BLCK).MODEL(M)=2.1;
                    continue
                case 'DRUCKER_PRAGER_I'
                    MATERIAL(BLCK).MODEL(M)=2.2;
                    continue
                case 'DRUCKER_PRAGER_PS'
                    MATERIAL(BLCK).MODEL(M)=2.3;
                    continue
                case 'MODIFIED_CAM_CLAY'
                    MATERIAL(BLCK).MODEL(M)=3.0;
                    continue
                case 'MODIFIED_CAM_CLAY_VISCO'
                    MATERIAL(BLCK).MODEL(M)=3.1;
                    continue
                case 'PZ_FORWARD'
                    MATERIAL(BLCK).MODEL(M)=4.2;
                    continue
                case 'PZ_MODIFIEDEULER'
                    MATERIAL(BLCK).MODEL(M)=4.3;
                    continue
                case 'PZ_BACKWARD'
                    MATERIAL(BLCK).MODEL(M)=4.1;
                    continue
                case 'VISCO_DEGRADATION'
                    MATERIAL(BLCK).MODEL(M)=5.0;
                    MATERIAL(BLCK).MODEL(M,2)=3;
                    continue
                otherwise
                    disp('Error, no such material model!')
                    stop
            end
        end
%         if strcmp(s1,'X_RANGE')
%             d = str2double(bb{t});
%             if isnan(d)
%                 if strcmp(bb{t},'FULL')
%                     RANGE(1,M*2-1)=L;
%                 elseif strcmp(bb{t},'INI')
%                     RANGE(1,M*2-1)=L0;
%                 else
%                     disp('Error, wrong material range!')
%                     stop
%                 end
%             else
%                 RANGE(1,M*2-1)=d;
%             end
%             d = str2double(c{t});
%             if isnan(d)
%                 if strcmp(c{t},'FULL')
%                     RANGE(1,M*2)=L;
%                 elseif strcmp(c{t},'INI')
%                     RANGE(1,M*2)=L0;
%                 else
%                     disp('Error, wrong material range!')
%                     stop
%                 end
%             else
%                 RANGE(1,M*2)=d;
%             end
%             continue
%         end
%         if strcmp(s1,'Y_RANGE')
%             d = str2double(bb{t});
%             if isnan(d)
%                 if strcmp(bb{t},'FULL')
%                     RANGE(2,M*2-1)=H;
%                 elseif strcmp(bb{t},'INI')
%                     RANGE(2,M*2-1)=H0;
%                 else
%                     disp('Error, wrong material range!')
%                     stop
%                 end
%             else
%                 RANGE(2,M*2-1)=d;
%             end
%             d = str2double(c{t});
%             if isnan(d)
%                 if strcmp(c{t},'FULL')
%                     RANGE(2,M*2)=H;
%                 elseif strcmp(c{t},'INI')
%                     RANGE(2,M*2)=H0;
%                 else
%                     disp('Error, wrong material range!')
%                     stop
%                 end
%             else
%                 RANGE(2,M*2)=d;
%             end
%             continue
%         end
        if strcmp(s1, '//')
        	continue
        end
        switch s1
            case 'BODY'
                BODIES(M)=str2double(bb{t});
                continue
            case 'EIGENEROSION'
                MATERIAL(BLCK).MODEL(M,2)=1;
                continue
            case 'EIGENSOFTENING'
                MATERIAL(BLCK).MODEL(M,2)=2;
                continue
            case 'YOUNG'
                MATERIAL(BLCK).MAT(1,M)={bb2(t)};
                continue
            case 'POISSON'
                MATERIAL(BLCK).MAT(2,M)={bb2(t)};
                continue
            case 'DENSITY'
                MATERIAL(BLCK).MAT(3,M)={bb2(t)};
                continue
            case 'SHEAR_MODULUS'
                MATERIAL(BLCK).MAT(4,M)={bb2(t)};
                continue
            case 'GHAR'
                MATERIAL(BLCK).MAT(4,M)={bb2(t)};
                continue
            case 'LAME_CONSTANT'
                MATERIAL(BLCK).MAT(5,M)={bb2(t)};
                continue
            case 'BULK_MODULUS'
                MATERIAL(BLCK).MAT(29,M)={bb2(t)};
                continue
            case 'KHAR'
                MATERIAL(BLCK).MAT(29,M)={bb2(t)};
                continue
            case 'CONSTRAINED_MODULUS'
                MATERIAL(BLCK).MAT(17,M)={bb2(t)};
                continue
            case 'WAVE_SPEED'
                MATERIAL(BLCK).MAT(6,M)={bb2(t)};
                continue
            case 'YIELD_STRESS'
                MATERIAL(BLCK).MAT(7,M)={bb2(t)};
                continue
            case 'COHESION'
                MATERIAL(BLCK).MAT(7,M)={bb2(t)};
                continue
            case 'PRECONSOLIDATION'
                MATERIAL(BLCK).MAT(7,M)={bb2(t)};    % Preconsolidation, before OCR
                continue
            case 'INITIAL_PRESSURE'
                MATERIAL(BLCK).MAT(25,M)=bb(t);%str2double(bb{t});    % Initial pressure
                continue
            case 'P0'
                MATERIAL(BLCK).MAT(25,M)=bb(t);%str2double(bb{t});    % Initial pressure
                continue
            case 'HARDENING'
                MATERIAL(BLCK).MAT(8,M)={bb2(t)};
                continue
            case 'HARDENING_EXPONENT'
                MATERIAL(BLCK).MAT(9,M)={bb2(t)};
                continue
            case 'EPSILON0'
                MATERIAL(BLCK).MAT(10,M)={bb2(t)};
                continue
            case 'FRICTION_ANGLE'
                MATERIAL(BLCK).MAT(11,M)={bb2(t)* pi/180};
                continue
            case 'DILATANCY_ANGLE'
                MATERIAL(BLCK).MAT(12,M)={bb2(t)* pi/180};
                continue
            case 'VISCOSITY'
                MATERIAL(BLCK).MAT(13,M)={bb2(t)};
                continue
            case 'GAMMA0'
                MATERIAL(BLCK).MAT(13,M)={bb2(t)};
                continue
            case 'VISCOSITY_EXPONENT'
                MATERIAL(BLCK).MAT(14,M)={bb2(t)};
                continue
            case 'N'
                MATERIAL(BLCK).MAT(14,M)={bb2(t)};
                continue
            case 'PERMEABILITY'
                MATERIAL(BLCK).MAT(15,M)={bb2(t)};
                continue
            case 'POROSITY'
                MATERIAL(BLCK).MAT(16,M)={bb2(t)};
                continue
            case 'WATER_BULK_MODULUS'
                MATERIAL(BLCK).MAT(18,M)={bb2(t)};
                continue
            case 'CRITICAL_STATE_LINE'
                MATERIAL(BLCK).MAT(19,M)={bb2(t)};
                continue
            case 'MF'
                MATERIAL(BLCK).MAT(19,M)={bb2(t)};
                continue
            case 'ALPHA_PARAMETER'
                MATERIAL(BLCK).MAT(20,M)={bb2(t)};
                continue
            case 'LAMBDA'
                MATERIAL(BLCK).MAT(21,M)={bb2(t)};
                continue
            case 'KAPPA'
                MATERIAL(BLCK).MAT(22,M)={bb2(t)};
                continue
            case 'INITIAL_VOLUMETRIC_STRAIN'
                MATERIAL(BLCK).MAT(23,M)=bb(t);%str2double(bb{t});
                continue
            case 'INITIAL_DEVIATORIC_STRAIN'
                MATERIAL(BLCK).MAT(26,M)=bb(t);%str2double(bb{t});
                continue
            case 'OCR'
                MATERIAL(BLCK).MAT(24,M)={bb2(t)};
                continue
            case 'KS'
                ks=bb2(t);
                continue
            case 'KW'
                kw=bb2(t);
                continue
            case 'CREEP_INDEX'
                MATERIAL(BLCK).MAT(30,M)={bb2(t)};
                continue
            case 'REFERENCE_TIME'
                MATERIAL(BLCK).MAT(31,M)={bb2(t)};
                continue
            case 'MG'
                MATERIAL(BLCK).MAT(32,M)={bb2(t)};
                continue
            case 'ALPHA_F'
                MATERIAL(BLCK).MAT(33,M)={bb2(t)};
                continue
            case 'ALPHA_G'
                MATERIAL(BLCK).MAT(34,M)={bb2(t)};
                continue
            case 'BETA0'
                MATERIAL(BLCK).MAT(35,M)={bb2(t)};
                continue
            case 'BETA1'
                MATERIAL(BLCK).MAT(36,M)={bb2(t)};
                continue
            case 'H0'
                MATERIAL(BLCK).MAT(37,M)={bb2(t)};
                continue
            case 'GAMMA_HDM'
                MATERIAL(BLCK).MAT(38,M)={bb2(t)};
                continue
            case 'HU0'
                MATERIAL(BLCK).MAT(39,M)={bb2(t)};
                continue
            case 'GAMMA_U'
                MATERIAL(BLCK).MAT(40,M)={bb2(t)};
                continue
            case 'GAMMA_VOL'
                MATERIAL(BLCK).MAT(41,M)={bb2(t)};
                continue
            case 'WATER_DENSITY'
                MATERIAL(BLCK).MAT(42,M)={bb2(t)};
                continue
            case 'CEPS'
                MATERIAL(BLCK).MAT(43,M)={bb2(t)};
                continue
            case 'GC'
                MATERIAL(BLCK).MAT(44,M)={bb2(t)};
                continue
            case 'WC'
                MATERIAL(BLCK).MAT(45,M)={bb2(t)};
                continue
            case 'FT'
                MATERIAL(BLCK).MAT(46,M)={bb2(t)};
                continue
            case 'WC_P'
                MATERIAL(BLCK).MAT(47,M)={bb2(t)};
                continue
            case 'FT_P'
                MATERIAL(BLCK).MAT(48,M)={bb2(t)};
                continue
            case 'D'
                MATERIAL(BLCK).MAT(49,M)={bb2(t)};
                continue
            case 'LAMBDA0'
                MATERIAL(BLCK).MAT(50,M)={bb2(t)};
                continue
            case 'XI_VG'
                MATERIAL(BLCK).MAT(50,M)={bb2(t)};
                continue
            case 'LAMBDA1'
                MATERIAL(BLCK).MAT(51,M)={bb2(t)};
                continue
            case 'BETA_RW'
                MATERIAL(BLCK).MAT(52,M)={bb2(t)};
                continue
            case 'ALPHA_RW'
                MATERIAL(BLCK).MAT(53,M)={bb2(t)};
                continue
            case 'LAMBDAD'
                MATERIAL(BLCK).MAT(54,M)={bb2(t)};
                continue
            case 'ALPHA_VG'
                MATERIAL(BLCK).MAT(54,M)={bb2(t)};
                continue
            case 'XRD'
                MATERIAL(BLCK).MAT(55,M)={bb2(t)};
                continue 
            case 'P0_VG'
                MATERIAL(BLCK).MAT(55,M)={bb2(t)};
                continue 
            case 'YR'
                MATERIAL(BLCK).MAT(56,M)={bb2(t)};
                continue  
            case 'SWR'
                MATERIAL(BLCK).MAT(56,M)={bb2(t)};
                continue 
            case 'XRW'
                MATERIAL(BLCK).MAT(57,M)={bb2(t)};
                continue
            case 'BETAD'
                MATERIAL(BLCK).MAT(58,M)={bb2(t)};
                continue   
            case 'M_VG'
                MATERIAL(BLCK).MAT(58,M)={bb2(t)};
                continue             
            case 'BETAW'
                MATERIAL(BLCK).MAT(59,M)={bb2(t)};
                continue            
            case 'BETA1_SW'
                MATERIAL(BLCK).MAT(60,M)={bb2(t)};
                continue
            case 'N_VG'
                MATERIAL(BLCK).MAT(60,M)={bb2(t)};
                continue
            case 'RETENTION_CURVE'
                MATERIAL(BLCK).MAT(61,M)=bb(t);
                continue
            case 'TAU95'
                MATERIAL(BLCK).MAT(62,M)={bb2(t)};
                continue            
            case 'DELTA95'
                MATERIAL(BLCK).MAT(63,M)={bb2(t)};
                continue
            case 'XI95'
                MATERIAL(BLCK).MAT(64,M)={bb2(t)};
                continue
            otherwise
                fprintf('Error, no such material property: %s !! \n',s1)
                stop
        end

    end
    
    fclose(fid); 
    
    %localization(RANGE);
    
    bds=max(BODIES);
    SOLVER.BODIES=max(1,bds);
    
    for i=1:mats
        for j=1:GEOMETRY.mat_points
            if mat_e(j)==i
                GEOMETRY.body(j)=BODIES(i);
            elseif bds==0
                GEOMETRY.body(j)=1;
            end
        end
        
        if SOLVER.UW
            n=MATERIAL(BLCK).MAT{16,i};
            if n==0
                fprintf('Error, no porosity !! \n')
            end
            if ks==0
                ks=1e40;       % Pa
                MATERIAL(BLCK).MAT(27,i)={ks};
            else
                MATERIAL(BLCK).MAT(27,i)={ks};
            end
            if isempty(MATERIAL(BLCK).MAT{18,i})
                if kw==0
                    fprintf('Error, no water bulk modulus !! \n')
                else
                    MATERIAL(BLCK).MAT(28,i)={kw};
                    MATERIAL(BLCK).MAT(18,i)={1/(n/kw+(1-n)/ks)};
                end
            else
                if kw==0
                    MATERIAL(BLCK).MAT(28,i)={1/(1/MATERIAL(BLCK).MAT{18,i})+(1-n)/n/ks};
                else
                    fprintf(...
                    'Error, two different values of the water bulk modulus !! \n')
                end
            end
        end
        if MATERIAL(BLCK).MODEL(i)<3 || (MATERIAL(BLCK).MODEL(i)<6 && MATERIAL(BLCK).MODEL(i)>=5)
            MATERIAL(BLCK).MAT(:,i)=elastic_tools(MATERIAL(BLCK).MAT(:,i));
        end
        if MATERIAL(BLCK).MODEL(i)<3 && MATERIAL(BLCK).MODEL(i)>=2
            MATERIAL(BLCK).MAT(:,i)=dp_tools(MATERIAL(BLCK).MAT(:,i),MATERIAL(BLCK).MODEL(i));
        elseif MATERIAL(BLCK).MODEL(i)<4 && MATERIAL(BLCK).MODEL(i)>=3
            MATERIAL(BLCK).MAT(:,i)=mcc_tools(MATERIAL(BLCK).MAT(:,i));
        elseif (MATERIAL(BLCK).MODEL(i)<5 && MATERIAL(BLCK).MODEL(i)>=4 )||...
            (MATERIAL(BLCK).MODEL(i)<6 && MATERIAL(BLCK).MODEL(i)>=5)
            MATERIAL(BLCK).MAT(:,i)=dp_tools(MATERIAL(BLCK).MAT(:,i),MATERIAL(BLCK).MODEL(i));
        end
        
        % Non local failure
        if MATERIAL(BLCK).MODEL(i,2)==1
            if SOLVER.FRAC==1 || SOLVER.FRAC==0
                SOLVER.FRAC=1;
            else
                error('Two different fracture criteria')
            end
        elseif MATERIAL(BLCK).MODEL(i,2)==2
            if SOLVER.FRAC==2 || SOLVER.FRAC==0
                SOLVER.FRAC=2;
            else
                error('Two different fracture criteria')
            end
            
            if isempty(MATERIAL(BLCK).MAT{47,i})
                MATERIAL(BLCK).MAT(47,i)={0};
            end
        elseif MATERIAL(BLCK).MODEL(i,2)==3
            if SOLVER.FRAC==3 || SOLVER.FRAC==0
                SOLVER.FRAC=3;
            else
                error('Two different fracture criteria')
            end
            
            if isempty(MATERIAL(BLCK).MAT{47,i})
                MATERIAL(BLCK).MAT(47,i)={0};
            end
        end
        
    end

end

function Mat=elastic_tools(Mat)

    global SOLVER
    
    if isempty(Mat{4}) && isempty(Mat{5}) && isempty(Mat{29})
        Mat(4)={Mat{1}/2/(1+Mat{2})};             % G 
        Mat(5)={2*Mat{4}*Mat{2}/(1-2*Mat{2})};    % Lambda
        Mat(29)={Mat{5}+2*Mat{4}};                % K 
    elseif isempty(Mat{1}) && isempty(Mat{5}) && isempty(Mat{29})
        Mat(1)={2*Mat{4}*(1+Mat{2})};             % E
        Mat(5)={2*Mat{4}*Mat{2}/(1-2*Mat{2})};    % Lambda
        Mat(29)={Mat{5}+2*Mat{4}};                % K 
    elseif isempty(Mat{1}) && isempty(Mat{4}) && isempty(Mat{29})
        Mat(1)={Mat{5}*(1+Mat{2})*(1-2*Mat{2})/Mat{2}};% E
        Mat(4)={Mat{1}/2/(1+Mat{2})};             % G 
        Mat(29)={Mat{5}+2*Mat{4}};                % K
    elseif isempty(Mat{1}) && isempty(Mat{2}) && isempty(Mat{29})
        Mat(1)={Mat{4}*(3*Mat{5}+2*Mat{4})/(Mat{4}+Mat{5})}; % E
        Mat(2)={Mat{5}/2/(Mat{5}+Mat{4})}; %nu
        Mat(29)={Mat{5}+2*Mat{4}};                % K
    elseif isempty(Mat{1}) && isempty(Mat{5}) && isempty(Mat{2})
        Mat(5)={Mat{29}-2/3*Mat{4}};                         % Lambda
        Mat(1)={Mat{4}*(3*Mat{5}+2*Mat{4})/(Mat{4}+Mat{5})}; % E
        Mat(2)={Mat{5}/2/(Mat{5}+Mat{4})}; %nu
    elseif isempty(Mat{1}) && isempty(Mat{5}) && isempty(Mat{4})
        Mat(1)={3*Mat{29}*(1-2*Mat{2})};  % E
        Mat(4)={Mat{1}/2/(1+Mat{2})};   % G  
        Mat(5)={Mat{29}-2/3*Mat{4}};      % Lambda
    elseif isempty(Mat{1}) && isempty(Mat{2}) && isempty(Mat{4})
        Mat(4)={3/2*(Mat{29}-Mat{5})};     % G
        Mat(2)={Mat{5}/2/(Mat{5}+Mat{4})}; %nu
        Mat(1)={3*Mat{29}*(1-2*Mat{2})};  % E
    end
    
    
    % Wave velocity
    Mat(17)={Mat{1}*(1-Mat{2})/((1+Mat{2})*(1-2*Mat{2}))};
    if SOLVER.UW
        M=Mat{17}+Mat{18};
    else
        M=Mat{17};
    end
    Mat(6)={sqrt(M/Mat{3})};     % cp
    
end

function Mat=dp_tools(Mat,MODEL)

    if isempty(Mat{10})
        if isempty(Mat{8})
            Mat(10)={1e10};    
        elseif Mat{9}==1
            Mat(10)= {-Mat{7}/Mat{8}};
        else
            error('wrond hardening parameters');
        end
    end

    if isempty(Mat{8})
        Mat{8}=1;
        Mat{9}=1e10;
    else
        if isempty(Mat{9})
            Mat{9}=1;
        end
    end
    
    if isempty(Mat{19}) & isempty(Mat{11}) & ~(MODEL==2.0 || MODEL==5.0)
        disp('Error, no critical state line!')
        stop
    elseif isempty(Mat{19})
        Mat(19)={6*sin(Mat{11})/(3-sin(Mat{11}))};
    elseif isempty(Mat{11})
        Mat(11)={asin((3*Mat{19})/(6+Mat{19}))};
    end


end

function Mat=pz_tools(Mat)

    global SOLVER

    % Mf
    if isempty(Mat{19}) && isempty(Mat{11})
        disp('Error, no critical state line!')
        stop
    elseif isempty(Mat{19})
        Mat(19)={6*sin(Mat{11})/(3-sin(Mat{11}))};
    elseif isempty(Mat{11})
        Mat(11)={asin((3*Mat{19})/(6+Mat{19}))};
    end
    
    %Mg
    if isempty(Mat{32}) && isempty(Mat{12})
        disp('Error, no critical state line!')
        stop
    elseif isempty(Mat(32))
        Mat(32)={6*sin(Mat(12))/(3-sin(Mat(12)))};
    elseif isempty(Mat(12))
        Mat(12)={asin((3*Mat(32))/(6+Mat(32)))};
    end
        
    %H0
    if isempty(Mat{37})
        if isempty(Mat{21}) || isempty(Mat{22})
            disp('We need more PZ parameters!')
            stop
        else
            Mat(37)={1/(Mat{21}-Mat{22})};
        end
    end
    
    %Hu0
    if isempty(Mat{39})
        if isempty( Mat{37}) && isempty(Mat{22})
            disp('We need more PZ parameters!')
            stop
        else
            if isempty(Mat{22})
                Mat(39)=Mat(37);
            else
                Mat(39)={1/Mat{22}};
            end
        end
    end
    
%     %Khar Ghar
%     Mat(29)={-Mat(29}/Mat(25}};
%     Mat(4) ={-Mat(4}/Mat(25}};
    
    % Non zero
    if isempty(Mat{33})
        disp('We need more PZ parameters!')
        stop
    elseif isempty(Mat{34})
        Mat(34)=Mat(33); % Alpha g = alpha f
    end
    
    %%% E
    Mat(2)={(3*Mat{29}-2*Mat{4})/(6*Mat{29}+2*Mat{4})};
    Mat(1)={2*Mat{4}*(1+Mat{2})};
    % Wave velocity
    Mat(17)={Mat{1}*(1-Mat{2})/((1+Mat{2})*(1-2*Mat{2}))};
    if SOLVER.UW
        M=Mat{17}+Mat{18};
    else
        M=Mat{17};
    end
    Mat(6)={sqrt(M/Mat{3})};     % cp
    
    if isempty(Mat{41})
        Mat(41)={0};     % gamma vol
    end
    
end

function Mat=mcc_tools(Mat)

    global SOLVER

    % M
    if isempty(Mat{19}) & isempty(Mat{11})
        disp('Error, no critical state line!')
        stop
    elseif isempty(Mat{19})
        Mat(19)={6*sin(Mat{11})/(3-sin(Mat{11}))};
    elseif isempty(Mat{11})
        Mat(11)={asin((3*Mat{19})/(6+Mat{19}))};
    end
    
    % OCR
    if ~isnan(str2double(Mat{25})) & ~isnan(Mat{7})
        Mat(24) = {Mat{7}/str2double(Mat{25})};
    elseif ~isnan(str2double(Mat{25})) && ~isnan(Mat{24})
        Mat(7) = {str2double(Mat{25})*Mat{24}};
    elseif ~isnan(Mat{7}) & ~isnan(Mat{24})
        Mat(25) = {num2str(Mat{7}/Mat{24})};
    elseif ~isnan(Mat{7})
        Mat(25)={num2str(Mat{7})};
        Mat(24)={1};
    elseif ~isnan(str2double(Mat{25}))
        Mat(7)={str2double(Mat{25})};
        Mat(24)={1};
    elseif ~isnan(Mat{24})
        disp('We need more MCC parameters!')
        stop
    end
    
    if isempty(Mat{30})
        Mat(30)={0};
    end
    
    if isempty(Mat{31})
        Mat(31)={0};
    end
    
    if isempty(Mat{21}) | isempty(Mat{22})
        disp('We need more MCC parameters!')
        stop
    end
    
    %%% Elastic
    if isempty(Mat{29})
        Mat(29)={-str2double(Mat{25})/Mat{22}};
    end
    if isempty(Mat{1}) && isempty(Mat{5}) && isempty(Mat{2})
        Mat(5)={Mat{29}-2/3*Mat{4}};                         % Lambda
        Mat(1)={Mat{4}*(3*Mat{5}+2*Mat{4})/(Mat{4}+Mat{5})}; % E
        Mat(2)={Mat{5}/2/(Mat{5}+Mat{4})}; %nu
    end
    
    % Wave velocity
    Mat(17)={Mat{1}*(1-Mat{2})/((1+Mat{2})*(1-2*Mat{2}))};
    if SOLVER.UW
        M=Mat{17}+Mat{18};
    else
        M=Mat{17};
    end
    Mat(6)={sqrt(M/Mat{3})};     % cp
    
end

% function localization(R)
% 
%     global MATERIAL GEOMETRY
%     
%     for m=1:MATERIAL(BLCK).number
%         for nodo=1:GEOMETRY.nodes
%             tol=GEOMETRY.h_nds(nodo,1)/5;
%             if (GEOMETRY.x_0(nodo,2)>=R(2,m*2-1)-tol) ...
%                     && (GEOMETRY.x_0(nodo,2)<=R(2,m*2)+tol) 
%                 if (GEOMETRY.x_0(nodo,1)>=R(1,m*2-1)-tol) ...
%                         && (GEOMETRY.x_0(nodo,1)<=R(1,m*2)+tol)
%                     MATERIAL(BLCK).n(nodo)=m;
%                 end
%             end
%         end
% 
%         for e=1:GEOMETRY.mat_points
%             tol=GEOMETRY.h_ini(e,1)/5;
%             if (GEOMETRY.xg_0(e,2)>=R(2,m*2-1)-tol) ...
%                     && (GEOMETRY.xg_0(e,2)<=R(2,m*2)+tol)
%                 if (GEOMETRY.xg_0(e,1)>=R(1,m*2-1)-tol) ...
%                         && (GEOMETRY.xg_0(e,1)<=R(1,m*2)+tol)
%                     MATERIAL(BLCK).e(e)=m;
%                 end
%             end
%         end
%         
%     end
% 
% end
% 


