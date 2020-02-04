
function [Disp_field]=explicit_corrector...
    (STEP,MATRIX,Disp_field,load_t,int_Ft,gamma)

    global GEOMETRY SOLVER
    
    time_step=STEP.dt;
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    nodes=GEOMETRY.nodes;

    %% Allocate
    a0  = Disp_field.a;
    v0  = Disp_field.v;

    Inv=zeros(nodes*sp,nodes*sp);
    [vs,as,load_s,int_Fs]=deal(zeros(nodes*sp,2));
    if SOLVER.UW==1
        [vw,aw,load_w,int_Fw]=deal(zeros(nodes*sp,2));
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
    else
        int_Fs=int_Ft;
        load_s=load_t;
        as=a0;
        vs=v0;
    end
    
    
    %% 0. Boundary conditions
    [boundary,~,~]=calculate_boundaries(STEP);

    %%  1.  Solver   %%%%

    if SOLVER.UW==1
        % 5.1 W Solver
        residual=MATRIX.l_mass_w*(int_Fs(:,1)-int_Fs(:,2)) - ...
            (MATRIX.l_mass_w-MATRIX.l_mass)*(int_Fw(:,1)-int_Fw(:,2))...
            + time_step*MATRIX.l_mass*MATRIX.l_damp*aw(:,2)...
            + MATRIX.l_mass_w*(load_s(:,1)-load_s(:,2))...
            - MATRIX.l_mass*(load_w(:,1)-load_w(:,2));

        mass_mat=MATRIX.l_mass_w*MATRIX.l_mass_w-time_step*gamma*...
            MATRIX.l_mass*MATRIX.l_damp - MATRIX.l_mass*MATRIX.l_mass_wn;

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
        
    else
        residual=-(int_Fs(:,1)-int_Fs(:,2))+(load_s(:,1)-load_s(:,2));
    end

    % 5.2 U Solver

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

    as(:,1)=Inv*residual(:,1)+as(:,2);


    %% 2. Corrector
    for i=1:nodes
        for j=1:sp
            if (boundary((i-1)*df+j)==0)
                vs((i-1)*sp+j,1)=vs((i-1)*sp+j,1)+...
                    gamma*time_step*as((i-1)*sp+j,1);
            end
            if SOLVER.UW==1
                if (boundary((i-1)*df+sp+j)==0)
                    vw((i-1)*sp+j,1)=vw((i-1)*sp+j,1)+gamma*time_step*...
                    aw((i-1)*sp+j,1);
                end
            end
        end
    end
    
    %% 3. Assemble vectors
    if SOLVER.UW==1
        [Disp_field.v(:,1)]=join_vector(vs(:,1),vw(:,1));
        [Disp_field.a(:,1)]=join_vector(as(:,1),aw(:,1));
    else
        Disp_field.a=as;
        Disp_field.v=vs;
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

function [vt]=join_vector(v1,v2)

    global GEOMETRY
    
    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
        
    for j=1:GEOMETRY.nodes
        vt((j-1)*df+1:(j-1)*df+sp,1)=v1((j-1)*sp+1:j*sp,1);
        vt((j-1)*df+sp+1:j*df,1)=v2((j-1)*sp+1:j*sp,1);
    end

end