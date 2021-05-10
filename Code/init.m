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
    STEP.FAIL=0;
    
    
    %----------------------------------------------------------------------
    % DEFINITION OF LENGTHS
    %----------------------------------------------------------------------        
    l0          = GEOMETRY.df*GEOMETRY.nodes;
    l2          = GEOMETRY.s_dim;
    elements    = GEOMETRY.mat_points;
    nodes       = GEOMETRY.nodes;
    if SOLVER.SMALL==0
        l1          = GEOMETRY.f_dim;
    end

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
    Int_var.gamma  = zeros(elements,3);
    Int_var.epsv  = zeros(elements,3);
    Int_var.Sy  = zeros(elements,2);
    Int_var.Sy_r= zeros(elements,2);
    Int_var.H   = zeros(elements,2);
    Int_var.eta = zeros(elements,2);
    Int_var.P0  = zeros(elements,3);
       
    Mat_state.title='Material state: material point information';

    if SOLVER.SMALL==0
        Mat_state.F=zeros(l1*elements,2);
        Mat_state.Be=zeros(l1*elements,2);
    else
        Mat_state.Es=zeros(l2*elements,2);
        Mat_state.Es_e=zeros(l2*elements,2);
    end
    Mat_state.Sigma=zeros(l2*elements,3);
    Mat_state.fint=zeros(l0,2);
    
    if SOLVER.UW==1  || SOLVER.UW==4
        Mat_state.k=zeros(elements,1);
        Mat_state.pw=zeros(elements,3);
        if SOLVER.UW==4
            Mat_state.sw=ones(elements,2);
            Mat_state.cs=zeros(elements,2);
            Mat_state.Q=zeros(elements,2);
        end
        if SOLVER.SMALL==0
            Mat_state.Fw=zeros(l1*elements,2);
        else
            Mat_state.Esw=zeros(l2*elements,2);
        end
    elseif SOLVER.UW==2 || SOLVER.UW==3
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
    
    GLOBAL.pw(elements,dim)=0;
    GLOBAL.Ps(elements,dim)=0;
    GLOBAL.Qs(elements,dim)=0;
    
    if SOLVER.SMALL==0
        GLOBAL.F(l1*elements,dim) = 0;
        GLOBAL.Be(l1*elements,dim) = 0;
    end
    GLOBAL.Sigma(l2*elements,dim) = 0;
    GLOBAL.Es(l2*elements,dim) = 0;
    GLOBAL.Es_p(l2*elements,dim) = 0;
    GLOBAL.xg(elements*GEOMETRY.sp,dim)=0;
    if SOLVER.UW==1  || SOLVER.UW==4
        if SOLVER.SMALL==0
            GLOBAL.Fw(l1*elements,dim)=0;
        else
            GLOBAL.Esw(l2*elements,dim)=0;
        end
        if SOLVER.UW==4
            GLOBAL.sw(1:elements,STEP.ste_p+1:dim)=1;
            GLOBAL.Q(elements,dim)=0;
            GLOBAL.cs(elements,dim)=0;
        end
    elseif SOLVER.UW==2 || SOLVER.UW==3
        GLOBAL.dpw(GEOMETRY.sp*elements,dim)=0;
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
        % P0
        if isempty(mmat{25,mati})
            p0=0;
        else
            p0=str2double(mmat{25,mati});
        end
        ev0=str2double(mmat{23,mati});
        es0=str2double(mmat{26,mati});
        if ~isnan(p0) && (isempty(mmat{23,mati}) || isempty(mmat{26,mati}))
            if isempty(mmat{23,mati})
                ev0=0;
            end
            if isempty(mmat{26,mati})
                es0=0;
            end
            Int_var.P0(i,1:3)=[p0 ev0 es0];
        elseif isnan(p0) || isnan(es0) || isnan(ev0)
            Int_var.P0(i,1:3)=VECTORS.fill_p0(mmat{25,mati},GLOBAL,i,...
                STEP,mmat{26,mati},mmat{23,mati});                
        else
            Int_var.P0(i,1:3)=[p0 ev0 es0];
        end
        % Yield stress
        if MODEL(mati)>=2 && MODEL(mati)<5
            if isempty(mmat{7,mati})
               Int_var.Sy(i,2)=0;
            else
                Int_var.Sy(i,2) = -mmat{7,mati};
            end
            GLOBAL.Sy(i,1)=Int_var.Sy(i,2);
        end
        if SOLVER.UW>0
            Mat_state.k(i) = ...
                mmat{15,mati}/...
                mmat{42,mati}/VARIABLE.g;
            if SOLVER.UW==4
                krw=W_retention.K(i,STEP.BLCK,Mat_state.sw);
                Mat_state.k(i)=Mat_state.k(i)*krw;
            end
        end
    end
    
    %----------------------------------------------------------------------
    % FRACTURE
    %----------------------------------------------------------------------
    if SOLVER.FRAC
        Mat_state.w=zeros(elements,1);
        Mat_state.status=zeros(elements,2);
        GLOBAL.w(elements,dim) = 0;
        GLOBAL.status(elements,dim) = 0;
        if SOLVER.FRAC>1
            Mat_state.e_ini=zeros(elements,4);
            GLOBAL.e_ini(elements,dim+1) = 0;
        end
        
        [MAT_POINT]=FRAC.eps_nb(MAT_POINT,STEP.BLCK);
        
        GLOBAL.Force(dim,SOLVER.BODIES*GEOMETRY.sp)=0;
        GLOBAL.E.D(dim,1)=0;
        GLOBAL.E.W(dim,SOLVER.BODIES)=0;
        GLOBAL.E.K(dim,SOLVER.BODIES)=0;
        STEP.ENERGY.D=0;
    end
    
    %----------------------------------------------------------------------
    % VECTORS OF ZEROS and ONES or taken from FILE
    %----------------------------------------------------------------------
    
    if SOLVER.INIT_file~=0

        % MAT_POINT
        [phases,~]=size(SOLVER.PHASES);
        for i=1:phases
            [MAT_POINT{i}]=LIB.list2S(MAT_POINT{i},'J',GLOBAL.J(:,ste_p));
            [MAT_POINT{i}]=LIB.list2S(MAT_POINT{i},'xg',...
                reshape(GLOBAL.xg(:,ste_p),[elements,GEOMETRY.sp]));
        end
        
        [Disp_field,Mat_state,Int_var,stiff_mtx,MAT_POINT,STEP]=VECTORS.Update_ini(...
            STEP,GLOBAL,Disp_field,Mat_state,Int_var,MAT_POINT);
        
    else
        Disp_field.x_a=GEOMETRY.x_0;
        GLOBAL.xg(:,1)=reshape(GEOMETRY.xg_0,[GEOMETRY.sp*GEOMETRY.mat_points,1]);
        
        if SOLVER.SMALL==0
            [Mat_state]=ini_F(Mat_state,l1,elements,GEOMETRY.sp,SOLVER.UW);
        end
        
        GLOBAL.ste_p=1;
        GLOBAL.tp(1)=STEP.t;
        
        [~,GLOBAL.d(:,1),GLOBAL.v(:,1),~]=calculate_boundaries(STEP,0);
        
        Disp_field.d(:,1)=GLOBAL.d(:,1);
        Disp_field.v(:,1)=GLOBAL.v(:,1);
        Disp_field.a(:,1)=GLOBAL.a(:,1);
        
        [Mat_state,MAT_POINT]=update_strain(Disp_field.d,Mat_state,MAT_POINT,STEP);
        
        [GLOBAL,Mat_state,stiff_mtx,Int_var,MAT_POINT,STEP,Disp_field]=initial_constitutive...
        (GLOBAL,MAT_POINT,Mat_state,Int_var,STEP,Disp_field);
    
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
    
    if UW==1 || UW==4
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
    dim=dim_2+ste_p+BLOCKS;
    fprintf('%i plot steps\n',dim);
    fprintf('Save %i times before the final\n',round(dim/SOLVER.SAVE_F));
    
    STEP.BLCK=BLK;
    STEP.ste=ste_aux+stini;
    STEP.dt=dt;
    STEP.t=tb+tp0;
end

function [GLOBAL,Mat_state,stiff_mtx,Int_var,MAT_POINT,STEP,Disp_field]=...
    initial_constitutive(GLOBAL,MAT_POINT,Mat_state,Int_var,STEP,Disp_field)

    global GEOMETRY SOLVER  MATERIAL
        
    Kt=4;
    df=GEOMETRY.df;
    stiff_mtx = zeros(df*GEOMETRY.nodes);
    
    Mat=GEOMETRY.material;
    MAT=MATERIAL(STEP.BLCK).MAT;
    MODEL=MATERIAL(STEP.BLCK).MODEL;
    
    dims=GEOMETRY.s_dim;

    for e=1:GEOMETRY.mat_points
            
        if SOLVER.SMALL==0            
            P0      = Int_var.P0(e,1:3);

            % MATRIX B and F
            b1=exp(P0(2)/3)^2;
            Be=[b1 1e-15 0; 1e-15 b1 0; 0 0 b1];
            [be]=LIB.m2v(Be);
            for i=1:GEOMETRY.f_dim
                Mat_state.Be((e-1)*GEOMETRY.f_dim+i,2)=be(i,1);
            end 
        end

        if MODEL(Mat(e))>=4 && MODEL(Mat(e))<5
            Int_var.H(e,2)=MAT{37,Mat(e)};
        end
    end
    % Constitutive calculation
    [stiff_mtx,Int_var,Mat_state,STEP]=Constitutive.strain2const(...
            STEP,Mat_state,Int_var,stiff_mtx,Kt,MAT_POINT);
        
    % ----------------------------
    % Internal forces
    % ----------------------------
    [Mat_state]=Constitutive.internal_forces(MAT_POINT,Mat_state,STEP.BLCK);
    
    for e=1:GEOMETRY.mat_points        
        for i=1:dims
            %Mat_state.Sigma((e-1)*dims+i,3)=Mat_state.Sigma((e-1)*dims+i,1);
            GLOBAL.Sigma((e-1)*dims+i,1)=Mat_state.Sigma((e-1)*dims+i,1);
        end
        
        [GLOBAL.Ps(e,1),GLOBAL.Qs(e,1)]=...
            VECTORS.invar(Mat_state.Sigma(:,1),e);   %PRESSURE
       
        if SOLVER.UW>0
            Pw=SOLVER.INITIAL_PORE_PRESSURE;
            Mat_state.pw(e,3)=Pw;
            GLOBAL.pw(e,1)=Pw+Mat_state.pw(e,1);
        end
         
        % Store other vectors
        if MODEL(Mat(e))>=2
            Int_var.Sy(e,2)   =   Int_var.Sy(e,1);
            Int_var.Sy_r(e,2) =   Int_var.Sy_r(e,1);
            if MODEL(Mat(e))>=4
                Int_var.H(e,2) =   Int_var.H(e,1);
                Int_var.eta(e,2) = Int_var.eta(e,1);
                Int_var.epsv(e,2)= Int_var.epsv(e,1);
            end
        end
            
    end
    
    %STEP=Time_Scheme.step(STEP,Disp_field);
    [Disp_field,Mat_state,Int_var]=VECTORS.Update(...
                Disp_field,Mat_state,Int_var);
    


end

