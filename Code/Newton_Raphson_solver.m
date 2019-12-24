
function [Disp_field,Mat_state,MAT_POINT]=...
            Newton_Raphson_solver(ste,stiff_mtx,mass_mtx,damp_mtx,load_s,...
            MAT_POINT,Disp_field,Int_var,Mat_state,BLCK)
        
    global SOLVER GEOMETRY
    
    clear error du
        
    % Initial values
    x_a = Disp_field.x_a;
    d0  = Disp_field.d;
    a0  = Disp_field.a;
    v0  = Disp_field.v;
    du=zeros(GEOMETRY.df*GEOMETRY.nodes,SOLVER.NR_iterations(BLCK));
    
    [matrix]=Time_Scheme.matrix(mass_mtx,stiff_mtx,damp_mtx,ste,BLCK);
    [InvK,~]=apply_conditions(0,ste,matrix,0);

    % Solve Newton Raphson
    TOL=SOLVER.rel_tolerance(BLCK);
    
    if ste==SOLVER.step_ini(BLCK)
        TOL=TOL*10000;
    elseif ste<SOLVER.step_ini(BLCK)+3
        TOL=TOL*100;
    end
        
    NR1=SOLVER.NR(BLCK);
    error_nr=zeros(SOLVER.NR_iterations(BLCK),1);
    iter=0;
    a=1;
    while iter<SOLVER.NR_iterations(BLCK)
        
        iter=iter+1;

        % 1. Evaluate residual 
        [GT]=Time_Scheme.calculation(d0,a0,v0,Mat_state.fint,mass_mtx,...
            damp_mtx,load_s(:,1),load_s(:,2),ste,BLCK);

        [~,GT]=apply_conditions(min(iter,2),ste,matrix,GT);
        
        if isnan(GT)   %|| NORMErec(iter,1)>emax
            error('Fallo en el NR global \n');
            SOLVER.FAIL=1;
        else
            % 2. Check for convergence
            if iter==1
                GT0=norm(GT);
                TOL=max(SOLVER.abs_tolerance(BLCK)*GT0,TOL);
            else
                [CONVER,error_nr,a,iter]=LIB.convergence(GT,GT0,error_nr,...
                    TOL,iter,SOLVER.NR_iterations(BLCK),a);
                if CONVER==1     
                    break
                elseif a==1
                    if iter>SOLVER.NR_iterations(BLCK)/2 && NR1>1
                        NR1=max(1,floor(NR1/2));
                    end
                else
                    if norm(d0(:,1))<1e-12 || norm(du(:,iter))<1e-15
                        break;
                    end
                    
                    d0(:,1)=d0(:,1)-du(:,iter+1);
                    [Mat_state,MAT_POINT]=update_F(d0,Mat_state,MAT_POINT);
                    [~,Int_var,Mat_state]=...
                        Constitutive(0,ste,Int_var,Mat_state,MAT_POINT,BLCK);
                    [GT]=Time_Scheme.calculation(d0,a0,v0,Mat_state.fint,mass_mtx,...
                        damp_mtx,load_s(:,1),load_s(:,2),ste,BLCK);
                    [~,GT]=apply_conditions(min(iter,2),ste,matrix,GT);
                end
            end
        end

        % 3. Solve for displacement increment
        [du(:,iter+1)]=Time_Scheme.solver_1(InvK,GT,a0,v0,ste,BLCK);
        du(:,iter+1)=a*du(:,iter+1);
        d0(:,1)=d0(:,1)+du(:,iter+1);
        
            % 3.1 Check
            if isnan(du(:,iter+1))
                disp('Nan in displacements')
                SOLVER.FAIL=1;
            end
        
        
        
        for j=1:GEOMETRY.nodes
            for i=1:GEOMETRY.sp
                x_a(j,i)=GEOMETRY.x_0(j,i)+d0((j-1)*GEOMETRY.df+i,1);
            end
        end

        % 4. Update
        
        % 4.1 Deformation gradient
        [Mat_state,MAT_POINT]=update_F(d0,Mat_state,MAT_POINT);

        % 4.2 Shape function
        if SOLVER.REMAPPING
%         [B,near,p,gamma_,lam_LME,REMAP,wrap,EP]=LME_EP(jacobians,...
%             volume,x_a,xg,B,near,p,gamma_,lam_LME,wrap,EP,ste);
        end
        
        % 5. Evaluate tangent matrix for NR iteration
        if rem(iter,NR1)==0

             if ste==1 || a<1e-4
                 [stiff_mtx,Int_var,Mat_state]=...
                 Constitutive(4,ste,Int_var,Mat_state,MAT_POINT,BLCK);
             else
                [stiff_mtx,Int_var,Mat_state]=...
                Constitutive(1,ste,Int_var,Mat_state,MAT_POINT,BLCK);
            end
            [matrix]=Time_Scheme.matrix(mass_mtx,stiff_mtx,damp_mtx,ste,BLCK);
            [InvK,~]=apply_conditions(0,ste,matrix,0);
        else       
            [~,Int_var,Mat_state]=...
                Constitutive(0,ste,Int_var,Mat_state,MAT_POINT,BLCK);
        end


    end
    
    Disp_field.x_a  = x_a;
    Disp_field.d    = d0;
        
end
