classdef Time_Scheme
    properties      
        af=0;
        am=0;
        delta=0;
        alpha=0;
        theta=1;
        gamma=0;
        beta=0;
        %tp;
        %t;
    end
    methods
        function obj=Time_Scheme(TIS,af,am,delta,alpha,theta,rho,i)
         
            global SOLVER
                     
            if SOLVER.DYN(i)==0 && TIS~=0
                disp('error with the time integration scheme,')
                disp('changed to STATIC, beta and gamma zero');
                %TIS(i)=0;
                alpha=0;
                delta=0;
            elseif TIS==1 && alpha~=0
                disp('error time integration scheme, Newmark 1')
                disp('no inertial terms: beta changed to zero')
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
%                      if TIS==3                   %GENERALIZED ALPHA
%                         delta=0.5+af-am;
%                         alpha=0.25*(1-am+af)^2;
%                     elseif TIS==4                  %HHT
%                         delta=(1+2*af)/2;
%                         alpha=0.25*(1+af)^2;      
%                     elseif TIS==6                   %WBZ
%                         delta=0.5-am;
%                         alpha=0.25*(1-am)^2;
%                      end
                end
            elseif TIS==5                  %WILSON-THETA
                 if SOLVER.IMPLICIT(i)==0
                    disp('error with the explicit time integration scheme,')
                    stop
                 end
                if theta==0
                    disp('Error, theta cannot be zero for Wilson')
                    stop
                end
                alpha=1/4;
                delta=1/2;
            elseif TIS==7                   %COLLOCATION METHOD
                if SOLVER.IMPLICIT(i)==0
                    disp('error with the explicit time integration scheme,')
                    stop
                end
                if theta==0 
                    disp('Error, theta cannot be zero for Collocation')
                    stop
                end
                if alpha==0 
                    disp('Error, alpha cannot be zero for Collocation')
                    stop
                end
                if delta==0 
                    disp('Error, delta cannot be zero for Collocation')
                    stop
                end
            end
            
            obj.af=af;
            obj.am=am;
            obj.delta=delta;
            obj.gamma=delta;
            obj.beta=alpha;
            obj.alpha=alpha;
            obj.theta=theta;
        end
     
    end
    methods(Static)
        function variables(i)
            global MATERIAL SOLVER GEOMETRY

            % CFL CALCULATION
            %------------------------------------------------
            h=GEOMETRY.h_ini;
            tt(GEOMETRY.mat_points,1)=0;
            material=GEOMETRY.material;
            for e=1:GEOMETRY.mat_points
                if SOLVER.UW
                    tt(e,1)=min(h(e)/MATERIAL(i).MAT{6,material(e)},...
                    h(e)/sqrt(MATERIAL(i).MAT{28,material(e)}/...
                    MATERIAL(i).MAT{42,material(e)}));
                else
                    tt(e,1)=h(e)/MATERIAL(i).MAT{6,material(e)};
                end
            end
            TT=min(tt);
            CFL=SOLVER.time_step(i)/TT;
            fprintf('%f of CFL\n',CFL);

            dt=SOLVER.time_step(i);
            tf=SOLVER.time_factor(i);
            
            %ste=floor(SOLVER.Time_final(i)/tf/dt)+1;
            tb=0;
            ste=0;
            while tb<SOLVER.Time_final(i)
                ste=ste+1;
                tb=tb+dt*tf;
                dt=dt*tf;
            end
            
            
            if i==1
                SOLVER.step_final(i)=ste;
                SOLVER.Time_final(i)=tb;
            else
                SOLVER.step_final(i)=ste+SOLVER.step_final(i-1);
                SOLVER.Time_final(i)=tb+SOLVER.Time_final(i-1);
            end
            
            if SOLVER.SAVE_I==1
                SOLVER.dim=floor(SOLVER.step_final(i)/SOLVER.SAVE_I);
            else
                SOLVER.dim=floor(SOLVER.step_final(i)/SOLVER.SAVE_I)+1;
            end
            
            fprintf('%i plot steps when finish BLOCK %i\n',SOLVER.dim,i);
            fprintf('Save %i times when finish BLOCK %i\n',...
                round(SOLVER.dim/SOLVER.SAVE_F),i);
        end
             
        function [GT]=calculation(d1,a1,v1,Fint,Mm,Cm,load,load1,STEP)
            
            global TIME SOLVER
            
            if SOLVER.DYN(STEP.BLCK)==1

                af=TIME{STEP.BLCK}.af;
                am=TIME{STEP.BLCK}.am;
                delta=TIME{STEP.BLCK}.delta;
                alpha=TIME{STEP.BLCK}.alpha;
                theta=TIME{STEP.BLCK}.theta;
                
                delta_t=STEP.dt;

                if alpha==0
                    A=0;
                    B=0;
                    C=0;
                    D=1/delta/delta_t;
                    E=1-1/delta;
                    F=0;
                else
                    A=(1-am)/alpha/delta_t/delta_t/theta/theta;
                    B=(1-am)/alpha/delta_t/theta;
                    C=(1-am)/2/alpha-1;
                    D=(1-af)*delta/alpha/delta_t/theta;
                    E=1-(1-af)*delta/alpha;
                    F=(1-af)*(1-delta/2/alpha)*delta_t*theta;
                end

                G=theta*(1-af);
                H=(1-af);

                du=d1(:,1)-d1(:,2);

                GT= G*(load-load1)+load1 - (H*(Fint(:,1)-Fint(:,2))+Fint(:,2))...
                        -Mm*(A*du-B*v1(:,2)-C*a1(:,2))...
                        -Cm*(D*du+E*v1(:,2)+F*a1(:,2));
            else
                GT= load - Fint(:,1);
            end

        end
        
        function [matrix]=matrix(mass_mtx,stiff_mtx,damp_mtx,STEP)
            
            global TIME SOLVER
                        
            if SOLVER.DYN(STEP.BLCK)==0
                matrix = stiff_mtx;  
            else

                af=TIME{STEP.BLCK}.af;
                am=TIME{STEP.BLCK}.am;
                delta=TIME{STEP.BLCK}.delta;
                alpha=TIME{STEP.BLCK}.alpha;
                theta=TIME{STEP.BLCK}.theta;
                
                delta_t=STEP.dt;

                if alpha==0
                    A=0;
                    D=1/delta/delta_t;
                else
                    A=(1-am)/alpha/delta_t/delta_t/theta/theta;
                    D=(1-af)*delta/alpha/delta_t/theta;
                end

                matrix = (1-af)*stiff_mtx + A*mass_mtx + D*damp_mtx;
                
            end
        end
        
        function [incr_d]=solver_1(InvK,FT,a1,v1,STEP)
            
            global TIME SOLVER
 
            if TIME{STEP.BLCK}.theta==1 || SOLVER.DYN(STEP.BLCK)
                incr_d=InvK*FT;
            else
                
                alpha=TIME{STEP.BLCK}.alpha;
                theta=TIME{STEP.BLCK}.theta;
                
                delta_t=STEP.dt;
                
                A=1/alpha/delta_t/delta_t/theta/theta;
                B=1/alpha/delta_t/theta;
                C=1/2/alpha;

                incr_d_th=InvK*FT;
                incr_a_th=A*incr_d_th-B*v1(:,2)-C*a1(:,2);

                incr_d=incr_d_th-(delta_t*(theta-1)*v1(:,2)+delta_t^2/2*(theta^2-1)*a1(:,2)...
                                  +alpha*delta_t^2*(theta^2-1)*incr_a_th);
            end

        %     a1(:,1)=a1(:,2)+incr_a_th/theta;
        %     v1(:,1)=v1(:,2)+delta_t*a1(:,2)+delta*delta_t/theta*incr_a_th;
        %     d1(:,1)=d1(:,2)+delta_t*v1(:,2)+delta_t^2/2*a1(:,2)+...
        %         alpha*delta_t^2/theta*incr_a_th;

        end
        
        function [a1,v1]=solver_2(d1,a1,v1,STEP)
            
            global TIME

            af=TIME{STEP.BLCK}.af;
            am=TIME{STEP.BLCK}.am;
            delta=TIME{STEP.BLCK}.delta;
            alpha=TIME{STEP.BLCK}.alpha;
            theta=TIME{STEP.BLCK}.theta;
            
            time_step=STEP.dt;

            du=d1(:,1)-d1(:,2);

            if alpha==0
                A=0;
                B=0;
                C=0;
                D=1/delta/time_step;
                E=1-1/delta;
                F=0;
            else
                A=(1-am)/alpha/time_step/time_step/theta/theta;
                B=(1-am)/alpha/time_step/theta;
                C=(1-am)/2/alpha-1;
                D=(1-af)*delta/alpha/time_step/theta;
                E=1-(1-af)*delta/alpha;
                F=(1-af)*(1-delta/2/alpha)*time_step*theta;
            end

            a1(:,1)=A*du-B*v1(:,2)-C*a1(:,2);
            v1(:,1)=D*du+E*v1(:,2)+F*a1(:,2);

        end
        
        function STEP=step(STEP,Disp_field)
            global SOLVER GEOMETRY TIME
            
            BLCK=STEP.BLCK;
            
            if SOLVER.autoadapt~=0
                tol=SOLVER.autoadapt(BLCK);
                sp=GEOMETRY.sp;
                df=GEOMETRY.df;
                nodes=GEOMETRY.nodes;
                alpha=TIME{BLCK}.alpha;
                d0  = Disp_field.d;
                if alpha
                    ob  = Disp_field.a;
                else
                    ob  = Disp_field.v;
                end
                if SOLVER.UW
                    [us,os]=deal(zeros(nodes*sp,2));
                    [uw,ow]=deal(zeros(nodes*(df-sp),2));
                    [us(:,2),uw(:,2)]=Time_Scheme.split_vector(d0(:,2),us(:,2),uw(:,2));
                    [us(:,1),uw(:,1)]=Time_Scheme.split_vector(d0(:,1),us(:,1),uw(:,1));
                    [os(:,2),ow(:,2)]=Time_Scheme.split_vector(ob(:,2),os(:,2),ow(:,2));
                    [os(:,1),ow(:,1)]=Time_Scheme.split_vector(ob(:,1),os(:,1),ow(:,1));
                    w1=os(:,1)-os(:,2);
                    w2=ow(:,1)-ow(:,2);
                    if alpha
                        eu1=STEP.dt^2*w1(:);
                        eu2=STEP.dt^2*w2(:);
                    else
                        eu1=STEP.dt*w1(:);
                        eu2=STEP.dt*w2(:);
                    end
                    error=max(norm(eu1)/norm(us),norm(eu2)/norm(uw));
                else
                    us=d0;
                    os=ob;
                    w1=os(:,1)-os(:,2);
                    if alpha
                        eu1=STEP.dt^2*w1(:);
                    else
                        eu1=STEP.dt*w1(:);
                    end
                    error=norm(eu1)/norm(us);
                end
                
                %Autoadaptative time step
                if error<0.2*tol || error>2*tol %
                    if alpha
                        STEP.dt = STEP.dt*(tol/error)^(1/3);
                    else
                        STEP.dt = STEP.dt*(tol/error)^(1/2);
                    end
                else
                    if STEP.error > SOLVER.r_tolerance(BLCK)
                        STEP.dt = STEP.dt*...
                            (SOLVER.r_tolerance(BLCK)/STEP.error)^(1/3);
                    else
                        STEP.dt = STEP.dt*SOLVER.time_factor(BLCK);
                    end
                end
            else
                STEP.dt = STEP.dt*SOLVER.time_factor(BLCK);
            end
            STEP.t  = STEP.t + STEP.dt;
        end
        
        function [v1,v2]=split_vector(vt,v1,v2)

            global GEOMETRY SOLVER

            df=GEOMETRY.df;
            sp=GEOMETRY.sp;

            for j=1:GEOMETRY.nodes
                v1((j-1)*sp+1:j*sp,1)=vt((j-1)*df+1:(j-1)*df+sp,1);
                if SOLVER.UW==1
                    v2((j-1)*sp+1:j*sp,1)=vt((j-1)*df+sp+1:j*df,1);
                elseif SOLVER.UW==2
                    v2(j,1)=vt(j*df,1);
                end   
            end

        end
        
    end
end


% delta=beta1=gamma
% alpha=beta2=beta
