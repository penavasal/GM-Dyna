
function [Disp_field,Mat_state,Shape_function,FAIL]=...
            Newton_Raphson_solver(ste,GT,InvK,mass_mtx,damp_mtx,load_s,...
            Shape_function,Disp_field,Int_var,Mat_state,FAIL)
        
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
    
    while error(iter) > TOL

        % 2.1 Solver        
        [du(:,iter)]=G_solver_1(InvK,GT,a0,v0,ste);
        d0(:,1)=d0(:,1)+du(:,iter);
            
        for j=1:GEOMETRY.nodes
            for i=1:GEOMETRY.sp
                x_a(j,i)=GEOMETRY.x_0(j,i)+d0((j-1)*GEOMETRY.df+i,1);
            end
         end

        % 2.2 Deformation gradient
        [Mat_state,Shape_function]=update_F(d0,Mat_state,Shape_function);

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

                [stiff_mtx,Int_var,Mat_state,FAIL]=...
                    Constitutive(1,ste,Int_var,Mat_state,Shape_function,FAIL);

                [matrix]=G_matrix(mass_mtx,stiff_mtx,damp_mtx,ste);

                [InvK,~]=apply_conditions(0,ste,matrix,0);
            else           
                [~,Int_var,Mat_state,FAIL]=...
                    Constitutive(0,ste,Int_var,Mat_state,Shape_function,FAIL);
            end
        end

        % 2.4 G calculation
        [GT]=G_calculation(d0,a0,v0,Mat_state.fint,mass_mtx,...
            damp_mtx,load_s(:,1),load_s(:,2),ste);

        [~,GT]=apply_conditions(2,ste,0,GT);

        % 2.5 Error
        iter=iter+1;
        error(iter)=norm(GT);
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
