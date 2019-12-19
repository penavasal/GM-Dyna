function [BLK,ste,ste_p,MAT_POINT,Disp_field,Int_var,Mat_state,GLOBAL,...
        stiff_mtx]=init(MAT_POINT)
    
    global GEOMETRY VARIABLE MATERIAL TIME SOLVER
    
    
    % Which is the initial time and step
    
    if SOLVER.INIT_file~=0
        s=strcat('load("',SOLVER.INIT_file,'.mat","GLOBAL")');
        eval(s);
        
        if SOLVER.INIT_STEP~=0
            ste_p=SOLVER.INIT_STEP;
        else
            ste_p=GLOBAL.ste_p;
        end
        
        [dim,ste,BLK]=fix_time(GLOBAL.tp,ste_p);
        SOLVER.dim=dim;
 
    else
        dim = SOLVER.dim;
        ste=1;
        ste_p=1;
        BLK=1;
    end
    
    
    %----------------------------------------------------------------------
    % DEFINITION OF LENGTHS
    %----------------------------------------------------------------------        
    l0          = GEOMETRY.df*GEOMETRY.nodes;
    l1          = GEOMETRY.f_dim;
    l2          = GEOMETRY.s_dim;
    elements    = GEOMETRY.mat_points;
    nodes       = GEOMETRY.nodes;

    %----------------------------------------------------------------------
    % DEFINITION OF STRUCTURES
    %----------------------------------------------------------------------
    
    Disp_field.title='Displacement field, velocity, acceleration and position';
    Disp_field.d  = zeros(l0,2);
    Disp_field.v  = zeros(l0,2);
    Disp_field.a  = zeros(l0,2);
    Disp_field.x_a= zeros(nodes,2);
    
    Int_var.title='Internal variables: plastic information';
    Int_var.dgamma  = zeros(elements,2);
    Int_var.gamma  = zeros(elements,2);
    Int_var.epsv  = zeros(elements,2);
    Int_var.Sy  = zeros(elements,2);
    Int_var.Sy_r= zeros(elements,2);
    Int_var.H   = zeros(elements,2);
    Int_var.eta = zeros(elements,2);
    Int_var.P0  = zeros(elements,1);
       
    Mat_state.title='Material state: material point information';

    Mat_state.F=zeros(l1*elements,2);
    Mat_state.Be=zeros(l1*elements,2);
    Mat_state.Sigma=zeros(l2*elements,2);
    Mat_state.fint=zeros(l0,2);
    
    if SOLVER.UW==1
        Mat_state.k=zeros(elements,1);
        Mat_state.pw=zeros(elements,3);
        Mat_state.Fw=zeros(l1*elements,2);
    elseif SOLVER.UW==2
        Mat_state.k=zeros(elements,1);
        Mat_state.pw=zeros(elements,3);
        Mat_state.dpw=zeros(GEOMETRY.sp*elements,2);
    end
    
    %----------------------------------------------------------------------
    % GLOBAL VECTORS
    %----------------------------------------------------------------------
    GLOBAL.title='Variables to save';
    GLOBAL.d(l0,dim)=0;
    GLOBAL.v(l0,dim)=0;
    GLOBAL.a(l0,dim)=0;
    
    GLOBAL.OutputList(dim,length(SOLVER.OutputType))=0;
    
    GLOBAL.dgamma(elements,dim)=0;
    GLOBAL.gamma(elements,dim)=0;
    GLOBAL.epsv(elements,dim)=0;
    GLOBAL.gamma_nds    = zeros(nodes,dim);
    GLOBAL.Sy(elements,dim)=0;
    GLOBAL.Sy_r(elements,dim)=0;
    GLOBAL.H(elements,dim)=0;
    GLOBAL.eta(elements,dim)=0;
    
    GLOBAL.Ps(elements,dim)=0;
    GLOBAL.Qs(elements,dim)=0;
       
    GLOBAL.F(l1*elements,dim) = 0;
    GLOBAL.Be(l1*elements,dim) = 0;
    GLOBAL.Sigma(l2*elements,dim) = 0;
    GLOBAL.Es(l2*elements,dim) = 0;
    GLOBAL.Es_p(l2*elements,dim) = 0;
    GLOBAL.xg(elements*GEOMETRY.sp,dim)=0;
    if SOLVER.UW==1
        GLOBAL.pw(elements,dim)=0;
        GLOBAL.Fw(l1*elements,dim)=0;
    end
    
    GLOBAL.tp(dim,1)=0;
    GLOBAL.J(1:elements,ste_p:dim)=1;
    
    %----------------------------------------------------------------------
    % VECTORS OF PARAMETERS
    %----------------------------------------------------------------------
    mmat=MATERIAL(BLK).MAT;
    for i=1:elements
        mati=GEOMETRY.material(i);
        Int_var.Sy(i,2) = mmat(7,mati);
        if SOLVER.UW>0
            Mat_state.k(i) = ...
                mmat(15,mati)/...
                mmat(42,mati)/VARIABLE.g;
        end
    end
    
    %----------------------------------------------------------------------
    % VECTORS OF ZEROS and ONES or taken from FILE
    %----------------------------------------------------------------------
    
    if SOLVER.INIT_file~=0

        for j=1:GEOMETRY.nodes
            for i=1:GEOMETRY.sp
                 Disp_field.x_a(j,i)=GEOMETRY.x_0(j,i)+...
                     GLOBAL.d((j-1)*GEOMETRY.df+i,ste_p);
            end
         end
                   
        Disp_field.a(:,2)=GLOBAL.a(:,ste_p);
        Disp_field.v(:,2)=GLOBAL.v(:,ste_p);
        Disp_field.d(:,2)=GLOBAL.d(:,ste_p);
        Disp_field.d(:,1)=GLOBAL.d(:,ste_p);

        Int_var.gamma(:,2)=GLOBAL.gamma(:,ste_p);
        Int_var.epsv(:,2)=GLOBAL.epsv(:,ste_p);
        Int_var.Sy(:,2)=GLOBAL.Sy(:,ste_p);
        Int_var.Sy_r(:,2)=GLOBAL.Sy_r(:,ste_p);
        Int_var.eta(:,2)=GLOBAL.eta(:,ste_p);
        Int_var.H(:,2)=GLOBAL.H(:,ste_p);
        Int_var.dgamma(:,2)=GLOBAL.dgamma(:,ste_p);
        Int_var.P0(:,1)=GLOBAL.Ps(:,1);
        
        Mat_state.Sigma(:,2)=GLOBAL.Sigma(:,ste_p);                %STRESS
        Mat_state.Sigma(:,3)=GLOBAL.Sigma(:,1);            %STRESS INITIAL
        Mat_state.F(:,2)=GLOBAL.F(:,ste_p);          %DEFORMATION GRADIENT
        Mat_state.F(:,1)=GLOBAL.F(:,ste_p);          %DEFORMATION GRADIENT
        Mat_state.Be(:,2)=GLOBAL.Be(:,ste_p);           %LEFT CAUCHY GREEN 
        
        [MAT_POINT]=LIB.list2S(MAT_POINT,'J',GLOBAL.J(:,ste_p));
        
        [MAT_POINT]=LIB.list2S(MAT_POINT,'xg',...
            reshape(GLOBAL.xg(:,ste_p),[elements,GEOMETRY.sp]));
        
        if SOLVER.UW>0
            Mat_state.pw(:,2)=GLOBAL.pw(:,ste_p);
            Mat_state.pw(:,3)=GLOBAL.pw(:,1);
            if SOLVER.UW==1
                Mat_state.Fw(:,2)=GLOBAL.Fw(:,ste_p);
                Mat_state.Fw(:,1)=GLOBAL.Fw(:,ste_p);
            end
        end
        
        MAT_POINT=shape_function_calculation(0,MAT_POINT,Disp_field);
        
        [stiff_mtx,Int_var,Mat_state]=...
        Constitutive(1,ste,Int_var,Mat_state,MAT_POINT,BLK);
        
        Mat_state.fint(:,2)=Mat_state.fint(:,1);
        
    else
        Disp_field.x_a=GEOMETRY.x_0;
        GLOBAL.xg(:,1)=reshape(GEOMETRY.xg_0,[GEOMETRY.sp*GEOMETRY.mat_points,1]);
        
        [Mat_state]=ini_F(Mat_state,l1,elements,GEOMETRY.sp,SOLVER.UW);
        
        GLOBAL.ste_p=1;
        GLOBAL.tp(1)=TIME{1}.t(1);
        
        [~,GLOBAL.d(:,1),GLOBAL.v(:,1),~]=calculate_boundaries(1,0);
        
        Disp_field.v(:,2)=GLOBAL.v(:,1);
        Disp_field.a(:,2)=GLOBAL.a(:,1);
        Disp_field.d(:,2)=GLOBAL.d(:,1);
        
        [GLOBAL,Mat_state,stiff_mtx,Int_var,MAT_POINT]=initial_constitutive...
        (GLOBAL,MAT_POINT,Mat_state,Int_var,1);
    end 

end

function [Mat_state]=ini_F(Mat_state,l1,elements,sp,UW)

    for e=1:elements       
        Mat_state.F((e-1)*l1+1,1)=1; 
        Mat_state.F((e-1)*l1+sp+2,1)=1;
        Mat_state.F(e*l1,1)=1;
        
        Mat_state.Be((e-1)*l1+1,1)=1; 
        Mat_state.Be((e-1)*l1+sp+2,1)=1;
        Mat_state.Be(e*l1,1)=1;  
        
    end
    Mat_state.F(:,2)=Mat_state.F(:,1);
    Mat_state.Be(:,2)=Mat_state.Be(:,1);
    
    if UW==1
        for e=1:elements       
            Mat_state.Fw((e-1)*l1+1,1)=1; 
            Mat_state.Fw((e-1)*l1+sp+2,1)=1;
            Mat_state.Fw(e*l1,1,1)=1;
        end
        Mat_state.Fw(:,2)=Mat_state.Fw(:,1);
    end
end

function [dim,ste_aux,BLK]=fix_time(tp,ste_p)

    global SOLVER TIME
    
    BLOCKS=SOLVER.BLOCKS;
    
    t_p=tp(ste_p);
    
    BLK=0;
    for i=1:BLOCKS
        if t_p<SOLVER.Time_final(i)
            BLK=i;
        end
    end
    if BLK==0
        disp('Init File time larger than problem time')
        stop
    end
    
    l1=length(TIME{BLK}.t);

    [~,ste_aux]=min(abs(TIME{BLK}.t(:)-t_p*ones(l1,1)));
    if TIME{BLK}.t(ste_aux)-t_p < 0
        ste_aux=ste_aux+1;
    end
    l2=l1-ste_aux;
    
    dim_2=floor((l2)/SOLVER.SAVE_I);
    dim=dim_2+ste_p;
    fprintf('%i plot steps\n',SOLVER.dim);
    fprintf('Save %i times before the final\n',round(SOLVER.dim/SOLVER.SAVE_F));
end

function [GLOBAL,Mat_state,stiff_mtx,Int_var,MAT_POINT]=...
    initial_constitutive(GLOBAL,MAT_POINT,Mat_state,Int_var,BLCK)

    global GEOMETRY SOLVER  MATERIAL
    
    MODEL=MATERIAL(BLCK).MODEL;
    Mat=GEOMETRY.material;
    MAT=MATERIAL(BLCK).MAT;
    
    df=GEOMETRY.df;
    dimf=GEOMETRY.f_dim;
    dims=GEOMETRY.s_dim;
        
    Kt=4;

    stiff_mtx = zeros(df*GEOMETRY.nodes);
    
    f_v       = zeros(dimf,1);
    f_old     = zeros(dimf,1);
    
    
    for e=1:GEOMETRY.mat_points

        % MATRIX B and F
        b1=exp(MAT(23,Mat(e))/3)^2;
        Be=[b1 1e-15 0; 1e-15 b1 0; 0 0 b1];
        
        for i=1:dimf
            f_v(i,1)=Mat_state.F((e-1)*dimf + i,1);
            f_old(i,1)=Mat_state.F((e-1)*dimf + i,2);
        end           
        [F]=LIB.v2m(f_v);
        [Fold]=LIB.v2m(f_old);
        jacobians=det(F);
        
            
        if MODEL(Mat(e))>=2
            Sy=Int_var.Sy(e,2);
            Sy_r=Sy;
            Gamma=Int_var.gamma(e,2);
            dgamma=Int_var.dgamma(e,2); 
        end

        if MODEL(Mat(e))<2
            if MODEL(Mat(e))==0
                [A,T,Be]=Saint_Venant(Kt,e,F,BLCK);
            elseif MODEL(Mat(e))<2 && MODEL(Mat(e))>=1
                [A,T,Be]=Neo_Hookean(Kt,e,F,jacobians,BLCK);
            end
        else        
            if MODEL(Mat(e))>=2 && MODEL(Mat(e))<3
                [A,T,Gamma,dgamma,Sy,Be]=...
                    Drucker_prager(Kt,e,Gamma,dgamma,Sy,F,Be,Fold,BLCK);
            elseif MODEL(Mat(e))>=3 && MODEL(Mat(e))<4
                [A,T,Gamma,dgamma,Sy,Sy_r,Be]=...
                    M_Cam_Clay(Kt,1,e,Gamma,dgamma,Sy,Sy_r,F,Fold,Be,press,BLCK);
                Int_var.Sy_r(e,1) = Sy_r;
            elseif MODEL(Mat(e))>=4 && MODEL(Mat(e))<5
                H=MAT(37,Mat(e));
                epsvol=0;
                [A,T,Gamma,epsvol,dgamma,Sy,etaB,H,Be]=...
                    PZ(Kt,1,e,Gamma,epsvol,dgamma,Sy,0,H,F,Fold,Be,press,BLCK);
                Int_var.epsv(e,1)  = epsvol;
            end
            Int_var.Sy(e,1)     = Sy;
            Int_var.gamma(e,1)  = Gamma;
            Int_var.dgamma(e,1) = dgamma;
        end

        % Main vectors
        [T_vec]=LIB.E2e(T);
        for i=1:dims
            Mat_state.Sigma((e-1)*dims+i,1)=T_vec(i,1);
            Mat_state.Sigma((e-1)*dims+i,3)=T_vec(i,1);
            GLOBAL.Sigma((e-1)*dims+i,1)=T_vec(i,1);
        end
        [be]=LIB.m2v(Be);
        [f]=LIB.m2v(F);
        for i=1:dimf
            Mat_state.Be((e-1)*dimf+i,2)=be(i,1);
            Mat_state.F((e-1)*dimf+i,2)=f(i,1);
        end 
        MAT_POINT(e).J    =   jacobians;
        
        if SOLVER.UW>0
            Pw=SOLVER.INITIAL_PORE_PRESSURE;
            Mat_state.pw(e,1)=Pw;
            Mat_state.pw(e,3)=Pw;
            GLOBAL.pw(e,1)=Pw;
        end
         
        % Stiffness matrix
        [stiff_mtx]=stiff_mat(MAT_POINT,Mat_state,e,stiff_mtx,T,A,BLCK);

        % Store other vectors
        if MODEL(Mat(e))>=2
            Int_var.Sy(e,2)   =   Int_var.Sy(e,1);
            Int_var.Sy_r(e,2) =   Int_var.Sy_r(e,1);
            Int_var.P0(e,1)   =   press;
            if MODEL(Mat(e))>=4
                Int_var.H(e,1) =   H;
                Int_var.H(e,2) =   H;
                Int_var.eta(e,1) = etaB;
                Int_var.eta(e,2) = etaB;
                Int_var.epsv(e,2)= epsvol;
            end
        end
            
    end

end

