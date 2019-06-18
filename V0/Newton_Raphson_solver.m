
function [Disp_field,Mat_state,MAT_POINT]=...
            Newton_Raphson_solver(ste,GT,InvK,mass_mtx,damp_mtx,load_s,...
            MAT_POINT,Disp_field,Int_var,Mat_state)
        
    global SOLVER GEOMETRY
    
    clear error du
    
    TOL=max(SOLVER.rel_tolerance*norm(GT),SOLVER.abs_tolerance);
    
    du=zeros(GEOMETRY.df*GEOMETRY.nodes,SOLVER.NR_iterations);

    error=zeros(SOLVER.NR_iterations,1);
    
    error(1)=1e32;
    iter=1;
    
    x_a = Disp_field.x_a;
    d0  = Disp_field.d;
    a0  = Disp_field.a;
    v0  = Disp_field.v;

    NR1=SOLVER.NR;
    
    a=1;
    while error(iter) > TOL

        % 2.1 Solver
        IK_s=a*InvK;
        [du(:,iter)]=Time_Scheme.solver_1(IK_s,GT,a0,v0,ste);
        d0(:,1)=d0(:,1)+du(:,iter);
            
        for j=1:GEOMETRY.nodes
            for i=1:GEOMETRY.sp
                x_a(j,i)=GEOMETRY.x_0(j,i)+d0((j-1)*GEOMETRY.df+i,1);
            end
         end

        % 2.2 Deformation gradient
        [Mat_state,MAT_POINT]=update_F(d0,Mat_state,MAT_POINT);

        if SOLVER.REMAPPING
%         [B,near,p,gamma_,lam_LME,REMAP,wrap,EP]=LME_EP(jacobians,...
%             volume,x_a,xg,B,near,p,gamma_,lam_LME,wrap,EP,ste);
        end

        % 2.3. Internal forces &/o stiffness matrix
        if iter>SOLVER.NR_iterations
            fprintf('No convergence achieved in step %i \n',ste);
            break;
            %stop;
        else
            if rem(iter,NR1)==0

                if ste==1 || a<1e-6
                    if a<1e-6
                        a=1;
                    end
                    [stiff_mtx,Int_var,Mat_state]=...
                    Constitutive(4,ste,Int_var,Mat_state,MAT_POINT);
                else
                    [stiff_mtx,Int_var,Mat_state]=...
                    Constitutive(1,ste,Int_var,Mat_state,MAT_POINT);
                end
                [matrix]=Time_Scheme.matrix(mass_mtx,stiff_mtx,damp_mtx,ste);

                [InvK,~]=apply_conditions(0,ste,matrix,0);
            else           
                [~,Int_var,Mat_state]=...
                    Constitutive(0,ste,Int_var,Mat_state,MAT_POINT);
            end
        end

        % 2.4 G calculation
        [GT]=Time_Scheme.calculation(d0,a0,v0,Mat_state.fint,mass_mtx,...
            damp_mtx,load_s(:,1),load_s(:,2),ste);

        [~,GT]=apply_conditions(2,ste,0,GT);

        % 2.5 Error
        iter=iter+1;
        error(iter)=norm(GT);
        if iter>5
            if error(iter)>error(iter-1)
                f1=error(iter-1);
                f2=error(iter);
                a=a*a*f1/2/(f2+f1*a-f1);
                iter = iter-1;
                if a<1e-18
                	fprintf('No convergence achieved in step %i \n',ste);
                    stop
                end
            else
                a=1;
            end
        end
        if error(iter)>TOL && iter>10
            if iter>SOLVER.NR_iterations/2
                if std(error(iter-10:iter-1))<TOL*1000
                    break;
                end
                fprintf('iter %i \n',iter);
                if NR1>1
                    NR1=max(1,floor(NR1/2));
                end
            else
                if std(error(iter-10:iter-1))<TOL*10
                    break;
                end
            end
        end
    end
    
    Disp_field.x_a  = x_a;
    Disp_field.d    = d0;
        
end
