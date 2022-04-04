
function [Disp_field]=explicit_corrector...
    (STEP,MATRIX,Disp_field,load_t,int_Ft)

    global GEOMETRY SOLVER TIME
    
    time_step=STEP.dt;
    BLCK=STEP.BLCK;
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    nodes=GEOMETRY.nodes;
    
    gamma=TIME{BLCK}.gamma;
    beta=TIME{BLCK}.beta;
    af=TIME{BLCK}.af;
    am=TIME{BLCK}.am;

    %% Allocate
    a0  = Disp_field.a;
    v0  = Disp_field.v;
    d0  = Disp_field.d;


    Inv=zeros(nodes*sp,nodes*sp);
    [us,vs,as,load_s,int_Fs]=deal(zeros(nodes*sp,2));
    if SOLVER.UW==1
        [uw,vw,aw,load_w,int_Fw]=deal(zeros(nodes*sp,2));
        [us(:,2),uw(:,2)]=split_vector(d0(:,2),us(:,2),uw(:,2));
        [vs(:,2),vw(:,2)]=split_vector(v0(:,2),vs(:,2),vw(:,2));
        [vs(:,1),vw(:,1)]=split_vector(v0(:,1),vs(:,1),vw(:,1));
        [as(:,2),aw(:,2)]=split_vector(a0(:,2),as(:,2),aw(:,2));
        [load_s(:,2),load_w(:,2)]=...
            split_vector(load_t(:,2),load_s(:,2),load_w(:,2));
        [load_s(:,1),load_w(:,1)]=...
            split_vector(load_t(:,1),load_s(:,1),load_w(:,1)); 
        [int_Fs(:,2),int_Fw(:,2)]=...
            split_vector(int_Ft(:,2),int_Fs(:,2),int_Fw(:,2));
        [int_Fs(:,1),int_Fw(:,1)]=...
            split_vector(int_Ft(:,1),int_Fs(:,1),int_Fw(:,1)); 
    elseif SOLVER.UW==2
        [pw,vpw,apw,load_pw,int_Fw]=deal(zeros(nodes,2));
        [us(:,2),pw(:,2)]=split_vector_up(d0(:,2),us(:,2),pw(:,2));
        [us(:,1),pw(:,1)]=split_vector_up(d0(:,1),us(:,1),pw(:,1));
        [vs(:,2),~]=split_vector_up(v0(:,2),vs(:,2),vpw(:,2));
        [vs(:,1),vpw(:,1)]=split_vector_up(v0(:,1),vs(:,1),vpw(:,1));
        [as(:,2),~]=split_vector_up(a0(:,2),as(:,2),apw(:,2));
        [load_s(:,2),load_pw(:,2)]=...
            split_vector_up(load_t(:,2),load_s(:,2),load_pw(:,2));
        [load_s(:,1),load_pw(:,1)]=...
            split_vector_up(load_t(:,1),load_s(:,1),load_pw(:,1)); 
        [int_Fs(:,2),int_Fw(:,2)]=...
            split_vector_up(int_Ft(:,2),int_Fs(:,2),int_Fw(:,2));
        [int_Fs(:,1),int_Fw(:,1)]=...
            split_vector_up(int_Ft(:,1),int_Fs(:,1),int_Fw(:,1)); 
    elseif SOLVER.UW==3
        [uw,vw,aw,load_w,int_Fw]=deal(zeros(nodes*sp,2));
        [pw,vpw,apw,load_pw,int_Fp]=deal(zeros(nodes,2));
        [us(:,2),uw(:,2),pw(:,2)]=split_vector_uwp(d0(:,2),us(:,2),uw(:,2),pw(:,2));
        [us(:,1),uw(:,1),pw(:,1)]=split_vector_uwp(d0(:,1),us(:,1),uw(:,1),pw(:,1));
        [vs(:,2),vw(:,2),~]=split_vector_uwp(v0(:,2),vs(:,2),vw(:,2),vpw(:,2));
        [vs(:,1),vw(:,1),vpw(:,1)]=split_vector_uwp(v0(:,1),vs(:,1),vw(:,1),vpw(:,1));
        [as(:,2),aw(:,2),~]=split_vector_uwp(a0(:,2),as(:,2),aw(:,2),apw(:,2));
        [load_s(:,2),load_w(:,2),load_pw(:,2)]=...
            split_vector_uwp(load_t(:,2),load_s(:,2),load_w(:,2),load_pw(:,2));
        [load_s(:,1),load_w(:,1),load_pw(:,1)]=...
            split_vector_uwp(load_t(:,1),load_s(:,1),load_w(:,1),load_pw(:,1)); 
        [int_Fs(:,2),int_Fw(:,2),int_Fp(:,2)]=...
            split_vector_uwp(int_Ft(:,2),int_Fs(:,2),int_Fw(:,2),int_Fp(:,2));
        [int_Fs(:,1),int_Fw(:,1),int_Fp(:,1)]=...
            split_vector_uwp(int_Ft(:,1),int_Fs(:,1),int_Fw(:,1),int_Fp(:,1)); 
    else
        int_Fs=int_Ft;
        load_s=load_t;
        as=a0;
        us=d0;
        vs=v0;
    end
    
    
    %% 0. Boundary conditions
    [boundary,i_disp,velo]=calculate_boundaries(STEP);

    %%  1.  Solver   %%%%

    if SOLVER.UW==1
        % 1.1 W Solver
        residual=MATRIX.l_mass_w*(int_Fs(:,1)-int_Fs(:,2))  ...
            - (MATRIX.l_mass_w-MATRIX.l_mass)*(int_Fw(:,1)-int_Fw(:,2))...
            + time_step*MATRIX.l_mass*MATRIX.l_damp_w*aw(:,2)...
            + MATRIX.l_mass_w*(load_s(:,1)-load_s(:,2))...
            - MATRIX.l_mass*(load_w(:,1)-load_w(:,2));

        mass_mat=MATRIX.l_mass_w*MATRIX.l_mass_w-time_step*gamma*...
            MATRIX.l_mass*MATRIX.l_damp_w - MATRIX.l_mass*MATRIX.l_mass_wn;

        for i=1:nodes
            for j=1:sp
                if (boundary(i*df+1-j)~=0)
                    mass_mat(i*sp+1-j,i*sp+1-j)=1;
                    Inv(i*sp+1-j,i*sp+1-j)=1;
                    residual(i*sp+1-j,1)=...
                        (vw(i*sp+1-j,1)-vw(i*sp+1-j,2))/time_step;
                else
                    Inv(i*sp+1-j,i*sp+1-j)=1/mass_mat(i*sp+1-j,i*sp+1-j);
                end
            end
        end

        daw=Inv*residual(:,1);
        aw(:,1)=aw(:,2)+daw;

         residual=(int_Fs(:,1)-int_Fs(:,2))-(int_Fw(:,1)-int_Fw(:,2))-...
             MATRIX.l_mass_w*daw+(load_s(:,1)-load_s(:,2));

    elseif SOLVER.UW==3
        % 1.1 W Solver
        residual=-MATRIX.l_mass_w*(int_Fs(:,1)-int_Fs(:,2))  ...
            + MATRIX.l_mass*(int_Fw(:,1)-int_Fw(:,2))...
            - MATRIX.l_mass*MATRIX.l_damp_w*(vw(:,1)-vw(:,2))... %- MATRIX.l_mass*MATRIX.l_damp_w*vw(:,1)...
            + MATRIX.l_mass_w*(load_s(:,1)-load_s(:,2))...
            - MATRIX.l_mass*(load_w(:,1)-load_w(:,2));

        mass_mat = MATRIX.l_mass*MATRIX.l_mass_wn - ...
            MATRIX.l_mass_w*MATRIX.l_mass_w;

        for i=1:nodes
            for j=1:sp
                if (boundary(i*df-j)~=0)
                    mass_mat(i*sp+1-j,i*sp+1-j)=1;
                    Inv(i*sp+1-j,i*sp+1-j)=1;
                    residual(i*sp+1-j,1)=...
                        (vw(i*sp+1-j,1)-vw(i*sp+1-j,2))/time_step;
                else
                    Inv(i*sp+1-j,i*sp+1-j)=1/mass_mat(i*sp+1-j,i*sp+1-j);
                end
            end
        end

        daw=Inv*residual(:,1);
        aw(:,1)=aw(:,2)+daw;

         residual=(int_Fs(:,1)-int_Fs(:,2))-...
             MATRIX.l_mass_w*daw+(load_s(:,1)-load_s(:,2));

    elseif SOLVER.UW==2
        residual=(int_Fs(:,1)-int_Fs(:,2))+(load_s(:,1)-load_s(:,2))-MATRIX.l_damp*(vs(:,1)-vs(:,2));
    else
        residual=-(int_Fs(:,1)-int_Fs(:,2))+(load_s(:,1)-load_s(:,2))-MATRIX.l_damp*(vs(:,1)-vs(:,2));
    end

    % 1.2 U Solver

    for i=1:nodes
        for j=1:sp
            if (boundary((i-1)*df+j)~=0)
                Inv((i-1)*sp+j,(i-1)*sp+j)=1;
                residual((i-1)*sp+j,1)=(vs((i-1)*sp+j,1)-...
                    vs((i-1)*sp+j,2))/time_step-as((i-1)*sp+j,2);
            else
                Inv((i-1)*sp+j,(i-1)*sp+j)=...
                    1/MATRIX.l_mass((i-1)*sp+j,(i-1)*sp+j);
            end
        end
    end

    as(:,1)=1/(1-am)*Inv*residual(:,1)+as(:,2);


    %% 2. Corrector
    for i=1:nodes
        for j=1:sp
            if af==1
                if boundary((i-1)*df+j)==0
                    vs((i-1)*sp+j,1)=vs((i-1)*sp+j,2)+...
                        (1-gamma)*time_step*as((i-1)*sp+j,2);
                    us((i-1)*sp+j,1)=us((i-1)*sp+j,2)+time_step*vs((i-1)*sp+j,2)...
                        +(0.5-beta)*time_step*time_step*as((i-1)*sp+j,2);
                elseif boundary((i-1)*df+j)==1
                    us((i-1)*sp+j,1)=us((i-1)*sp+j,2)+i_disp((i-1)*sp+j,1);
                    vs((i-1)*sp+j,1)=i_disp((i-1)*sp+j,1)/time_step;
                else
                    vs((i-1)*sp+j,1)=velo((i-1)*sp+j,1);
                    us((i-1)*sp+j,1)=us((i-1)*sp+j,2)+velo((i-1)*sp+j,1)*time_step;
                end
            end
            
            
            if (boundary((i-1)*df+j)==0)
                vs((i-1)*sp+j,1)=vs((i-1)*sp+j,1)+...
                    gamma*time_step*as((i-1)*sp+j,1);
                us((i-1)*sp+j,1)=us((i-1)*sp+j,1)+...
                    beta*time_step*time_step*as((i-1)*sp+j,1);
            end
            if SOLVER.UW==1 || SOLVER.UW==3
                if (boundary((i-1)*df+sp+j)==0)
                    vw((i-1)*sp+j,1)=vw((i-1)*sp+j,1)+gamma*time_step*...
                    aw((i-1)*sp+j,1);
                    uw((i-1)*sp+j,1)=uw((i-1)*sp+j,1)+beta*time_step*...
                    time_step*aw((i-1)*sp+j,1);
                end
            end
        end
    end
    
    if SOLVER.UW==2 
        %2.1 pw_dot Corrector
        pw_dot = MATRIX.l_damp_w*vs(:,1) + load_pw(:,1) - int_Fw(:,1) + ...
            MATRIX.l_mass_w*as(:,1);
        for i=1:nodes
            if (boundary((i-1)*df+sp+1)==0)
                pw(i,1)=pw(i,1)+...
                    gamma*time_step*pw_dot(i);
                vpw(i,1)=pw_dot(i,1);
            end
        end
    elseif SOLVER.UW==3
        Inv=zeros(nodes);
        residual=int_Fp(:,1);
        for i=1:nodes
            if (boundary(i*df)~=0)
                Inv(i,i)=1;
                residual(i,1)=(pw(i,1)-...
                    pw(i,2))/time_step;%-vpw(i,2);
            else
                Inv(i,i)=...
                    1/MATRIX.l_mass_p(i,i);
            end
        end
        pw_dot=Inv*residual;
        for i=1:nodes
            if (boundary(i*df)==0)
                pw(i,1)=pw(i,1)+...
                    gamma*time_step*pw_dot(i);
                vpw(i,1)=pw_dot(i,1);
            end
        end
    end
    
    %% 3. Assemble vectors
    if SOLVER.UW==1
        [Disp_field.v(:,1)]=join_vector(vs(:,1),vw(:,1));
        [Disp_field.a(:,1)]=join_vector(as(:,1),aw(:,1));
    elseif SOLVER.UW==2
        [Disp_field.d(:,1)]=join_vector_up(us(:,1),pw(:,1));
        [Disp_field.v(:,1)]=join_vector_up(vs(:,1),vpw(:,1));
        [Disp_field.a(:,1)]=join_vector_up(as(:,1),zeros(nodes,1)); 
    elseif SOLVER.UW==3
        [Disp_field.d(:,1)]=join_vector_uwp(us(:,1),uw(:,1),pw(:,1));
        [Disp_field.v(:,1)]=join_vector_uwp(vs(:,1),vw(:,1),vpw(:,1));
        [Disp_field.a(:,1)]=join_vector_uwp(as(:,1),aw(:,1),zeros(nodes,1));
    else
        Disp_field.a=as;
        Disp_field.v=vs;
        Disp_field.d=us;
    end
    
end

function [v1,v2]=split_vector(vt,v1,v2)

    global GEOMETRY
    
    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
        
    for j=1:GEOMETRY.nodes
        v1((j-1)*sp+1:j*sp,1)=vt((j-1)*df+1:(j-1)*df+sp,1);
        v2((j-1)*sp+1:j*sp,1)=vt((j-1)*df+sp+1:j*df,1);
    end

end

function [v1,v2]=split_vector_up(vt,v1,v2)

    global GEOMETRY
    
    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
        
    for j=1:GEOMETRY.nodes
        v1((j-1)*sp+1:j*sp,1)=vt((j-1)*df+1:(j-1)*df+sp,1);
        v2(j,1)=vt(j*df,1);
    end

end

function [v1,v2,v3]=split_vector_uwp(vt,v1,v2,v3)

    global GEOMETRY
    
    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
        
    for j=1:GEOMETRY.nodes
        v1((j-1)*sp+1:j*sp,1)=vt((j-1)*df+1:(j-1)*df+sp,1);
        v2((j-1)*sp+1:j*sp,1)=vt((j-1)*df+sp+1:(j-1)*df+sp*2,1);
        v3(j,1)=vt(j*df,1);
    end

end

function [vt]=join_vector(v1,v2)

    global GEOMETRY
    
    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
        
    for j=1:GEOMETRY.nodes
        vt((j-1)*df+1:(j-1)*df+sp,1)=v1((j-1)*sp+1:j*sp,1);
        vt((j-1)*df+sp+1:j*df,1)=v2((j-1)*sp+1:j*sp,1);
    end

end

function [vt]=join_vector_up(v1,v2)

    global GEOMETRY
    
    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
    
    
    vt=zeros(df*GEOMETRY.nodes,1);
        
    for j=1:GEOMETRY.nodes
        vt((j-1)*df+1:(j-1)*df+sp,1)=v1((j-1)*sp+1:j*sp,1);
        vt(j*df,1)=v2(j,1);
    end

end

function [vt]=join_vector_uwp(v1,v2,v3)

    global GEOMETRY
    
    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
    
    
    vt=zeros(df*GEOMETRY.nodes,1);
        
    for j=1:GEOMETRY.nodes
        vt((j-1)*df+1:(j-1)*df+sp,1)=v1((j-1)*sp+1:j*sp,1);
        vt((j-1)*df+sp+1:j*df-1,1)=v2((j-1)*sp+1:j*sp,1);
        vt(j*df,1)=v3(j,1);
    end

end
