classdef Time_Scheme
    properties      
        af=0;
        am=0;
        delta=0;
        alpha=0;
        theta=1;
        gamma=0;
        tp;
        t;
    end
    methods
        function obj=Time_Scheme(TIS,af,am,delta,alpha,theta,rho)
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
            
            obj.af=af;
            obj.am=am;
            obj.delta=delta;
            obj.gamma=delta;
            obj.alpha=alpha;
            obj.theta=theta;
        end
     
    end
    methods(Static)
        function variables
            global MATERIAL VARIABLE SOLVER GEOMETRY TIME

            h=GEOMETRY.h_ini;
            tt(GEOMETRY.elements,1)=0;
            for e=1:GEOMETRY.elements
                if SOLVER.UW
                    tt(e,1)=min(h(e)/MATERIAL.MAT(6,MATERIAL.e(e)),...
                    h(e)/sqrt(MATERIAL.MAT(28,MATERIAL.e(e))/VARIABLE.rho_w));
                else
                    tt(e,1)=h(e)/MATERIAL.MAT(6,MATERIAL.e(e));
                end
            end
            TT=min(tt);
            CFL=SOLVER.time_step/TT;
            fprintf('%f of CFL\n',CFL);

            delta_t=SOLVER.time_step;

            ste=1;
            t(ste,1)=0;
            while t(ste,1)<SOLVER.Time_final
                ste=ste+1;
                delta_t=delta_t*SOLVER.time_factor;
                t(ste,1)=t(ste-1,1)+delta_t;
            end
            SOLVER.Time_final=t(ste,1);
            SOLVER.step_final=ste;
            TIME.t=t;


            if SOLVER.SAVE_I==1
                SOLVER.dim=floor(SOLVER.step_final/SOLVER.SAVE_I);
            else
                SOLVER.dim=floor(SOLVER.step_final/SOLVER.SAVE_I)+1;
            end
            TIME.tp=zeros(SOLVER.dim,1);
            fprintf('%i plot steps\n',SOLVER.dim);
            fprintf('Save %i times before the final\n',round(SOLVER.dim/SOLVER.SAVE_F));

        end
             
        function [GT]=calculation(d1,a1,v1,Fint,Mm,Cm,load,load1,ste)
            
            global TIME

            af=TIME.af;
            am=TIME.am;
            delta=TIME.delta;
            alpha=TIME.alpha;
            theta=TIME.theta;

            if ste==1
                delta_t=TIME.t(ste+1)-TIME.t(ste);
            else
                delta_t=TIME.t(ste)-TIME.t(ste-1);
            end

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

        end
        
        function [matrix]=matrix(mass_mtx,stiff_mtx,damp_mtx,ste)
            
            global TIME

            af=TIME.af;
            am=TIME.am;
            delta=TIME.delta;
            alpha=TIME.alpha;
            theta=TIME.theta;


            if ste==1
                delta_t=TIME.t(ste+1)-TIME.t(ste);
            else
                delta_t=TIME.t(ste)-TIME.t(ste-1);
            end

            if alpha==0
                A=0;
                D=1/delta/delta_t;
            else
                A=(1-am)/alpha/delta_t/delta_t/theta/theta;
                D=(1-af)*delta/alpha/delta_t/theta;
            end


            matrix = (1-af)*stiff_mtx + A*mass_mtx + D*damp_mtx;

        end
        
        function [incr_d]=solver_1(InvK,FT,a1,v1,ste)
            
            global TIME

            alpha=TIME.alpha;
            theta=TIME.theta;


            if ste==1
                delta_t=TIME.t(ste+1)-TIME.t(ste);
            else
                delta_t=TIME.t(ste)-TIME.t(ste-1);
            end

            if theta==1
                incr_d=InvK*FT;
            else
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
        
        function [a1,v1]=solver_2(d1,a1,v1,ste)
            
            global TIME

            af=TIME.af;
            am=TIME.am;
            delta=TIME.delta;
            alpha=TIME.alpha;
            theta=TIME.theta;


            du=d1(:,1)-d1(:,2);

            if ste==1
                time_step=TIME.t(ste+1)-TIME.t(ste);
            else
                time_step=TIME.t(ste)-TIME.t(ste-1);
            end

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
        
    end
end


% delta=beta1=gamma
% alpha=beta2=beta
