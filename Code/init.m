function [STEP,MAT_POINT,Disp_field,Int_var,Mat_state,GLOBAL,...
        stiff_mtx]=init(MAT_POINT)
    
    global GEOMETRY VARIABLE MATERIAL SOLVER
    
    
    % Which is the initial time and step
    
    if SOLVER.INIT_file~=0
        s=strcat('load("',SOLVER.INIT_file,'.mat","GLOBAL")');
        eval(s);
        
        if SOLVER.INIT_STEP~=0
            ste_p=SOLVER.INIT_STEP;
        else
            ste_p=GLOBAL.ste_p;
        end
        
        [dim,STEP]=fix_time(GLOBAL.tp,ste_p);
        SOLVER.dim=dim;
        STEP.ste_p=ste_p;
    else
        dim = SOLVER.dim;
        STEP.ste=1;
        STEP.ste_p=1;
        STEP.BLCK=1;
        STEP.t=0;
        STEP.dt=SOLVER.time_step(STEP.BLCK);
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
    Disp_field.d  = zeros(l0,3);
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
    GLOBAL.J(1:elements,STEP.ste_p+1:dim)=1;
    
    %----------------------------------------------------------------------
    % VECTORS OF PARAMETERS
    %----------------------------------------------------------------------
    mmat=MATERIAL(STEP.BLCK).MAT;
    MODEL=MATERIAL(STEP.BLCK).MODEL;

    for i=1:elements
        mati=GEOMETRY.material(i);
        if MODEL(mati)>=2 && MODEL(mati)<5
            % Yield stress
            if isempty(mmat{7,mati})
               Int_var.Sy(i,2)=0;
            else
                Int_var.Sy(i,2) = mmat{7,mati};
            end
            if MODEL(mati)>=3 && MODEL(mati)<5
                % P0
                p0=str2double(mmat{25,mati});
                ev0=str2double(mmat{23,mati});
                es0=str2double(mmat{26,mati});
                if isnan(p0)
                    Int_var.P0(i,1:3)=VECTORS.fill_p0(mmat{25,mati},GLOBAL,i,STEP);
                else
                    if isnan(ev0)
                        ev0=0;
                    end
                    if isnan(es0)
                        es0=0;
                    end                   
                    Int_var.P0(i,1:3)=[p0 ev0 es0]; 
                end
            end
        end
        if SOLVER.UW>0
            Mat_state.k(i) = ...
                mmat{15,mati}/...
                mmat{42,mati}/VARIABLE.g;
        end
    end
    
    %----------------------------------------------------------------------
    % VECTORS OF ZEROS and ONES or taken from FILE
    %----------------------------------------------------------------------
    
    if SOLVER.INIT_file~=0

        % MAT_POINT
        [MAT_POINT]=LIB.list2S(MAT_POINT,'J',GLOBAL.J(:,ste_p));
        [MAT_POINT]=LIB.list2S(MAT_POINT,'xg',...
            reshape(GLOBAL.xg(:,ste_p),[elements,GEOMETRY.sp]));
        MAT_POINT=shape_function_calculation(0,MAT_POINT,Disp_field);
        
        [Disp_field,Mat_state,Int_var,stiff_mtx]=VECTORS.Update_ini(...
            STEP,GLOBAL,Disp_field,Mat_state,Int_var,MAT_POINT);
        
    else
        Disp_field.x_a=GEOMETRY.x_0;
        GLOBAL.xg(:,1)=reshape(GEOMETRY.xg_0,[GEOMETRY.sp*GEOMETRY.mat_points,1]);
        
        [Mat_state]=ini_F(Mat_state,l1,elements,GEOMETRY.sp,SOLVER.UW);
        
        GLOBAL.ste_p=1;
        GLOBAL.tp(1)=STEP.t;
        
        [~,GLOBAL.d(:,1),GLOBAL.v(:,1),~]=calculate_boundaries(STEP,0);
        
        Disp_field.v(:,2)=GLOBAL.v(:,1);
        Disp_field.a(:,2)=GLOBAL.a(:,1);
        Disp_field.d(:,2)=GLOBAL.d(:,1);
        
        [GLOBAL,Mat_state,stiff_mtx,Int_var,MAT_POINT]=initial_constitutive...
        (GLOBAL,MAT_POINT,Mat_state,Int_var,STEP);
    
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

function [dim,STEP]=fix_time(tp,ste_p)

    global SOLVER
    
    BLOCKS=SOLVER.BLOCKS;
    
    t_p=tp(ste_p);
    
    BLK=0;
    for i=1:BLOCKS
        if t_p<SOLVER.Time_final(i)
            BLK=i;
            break;
        end
    end
    if BLK==0
        disp('Init File time larger than problem time')
        stop
    elseif BLK==1
        tp0=0;
        stini=0;
    else
        tp0=SOLVER.Time_final(BLK-1);
        stini=SOLVER.step_final(BLK-1);
    end
    
    dt=SOLVER.time_step(BLK);
    tf=SOLVER.time_factor(BLK);
    
    tb=0;
    ste=0;
    while tb<(t_p-tp0)
        ste=ste+1;
        tb=tb+dt*tf;
        dt=dt*tf;
    end
    tb=tb-dt;
    dt=dt/tf;
    ste_aux=ste-1;
    l1=SOLVER.step_final(BLOCKS);
    l2=l1-ste_aux;
    
    dim_2=floor((l2)/SOLVER.SAVE_I);
    dim=dim_2+ste_p;
    fprintf('%i plot steps\n',dim);
    fprintf('Save %i times before the final\n',round(dim/SOLVER.SAVE_F));
    
    STEP.BLCK=BLK;
    STEP.ste=ste_aux+stini;
    STEP.dt=dt;
    STEP.t=tb+tp0;
end

function [GLOBAL,Mat_state,stiff_mtx,Int_var,MAT_POINT]=...
    initial_constitutive(GLOBAL,MAT_POINT,Mat_state,Int_var,STEP)

    global GEOMETRY SOLVER  MATERIAL
    
    BLCK=STEP.BLCK;
    
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
        
        P0      = Int_var.P0(e,1:3);

        % MATRIX B and F
        b1=exp(P0(2)/3)^2;
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
                    M_Cam_Clay(Kt,1,e,Gamma,dgamma,Sy,Sy_r,F,Fold,Be,press,P0,BLCK);
                Int_var.Sy_r(e,1) = Sy_r;
            elseif MODEL(Mat(e))>=4 && MODEL(Mat(e))<5
                H=MAT(37,Mat(e));
                epsvol=0;
                [A,T,Gamma,epsvol,dgamma,Sy,etaB,H,Be]=...
                    PZ(Kt,1,e,Gamma,epsvol,dgamma,Sy,0,H,F,Fold,Be,P0,BLCK);
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
        
        [GLOBAL.Ps(e,1),GLOBAL.Qs(e,1)]=...
            VECTORS.invar(Mat_state.Sigma(:,1),e);   %PRESSURE
        
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
            %Int_var.P0(e,1)   =   press;
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

