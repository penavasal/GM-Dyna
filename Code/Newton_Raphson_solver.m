
function [Disp_field,Mat_state,MAT_POINT,STEP]=...
            Newton_Raphson_solver(STEP,stiff_mtx,mass_mtx,damp_mtx,load_s,...
            MAT_POINT,Disp_field,Int_var,Mat_state)
        
    global SOLVER GEOMETRY
    
    BLCK=STEP.BLCK;
    
    clear error du
        
    % Initial values
    x_a = Disp_field.x_a;
    d0  = Disp_field.d;
    a0  = Disp_field.a;
    v0  = Disp_field.v;
    du=zeros(GEOMETRY.df*GEOMETRY.nodes,SOLVER.NR_iterations(BLCK));
    
    [matrix]=Time_Scheme.matrix(mass_mtx,stiff_mtx,damp_mtx,STEP);
    [InvK,~]=apply_conditions(0,STEP,matrix,0);
    
    [GT0]=Time_Scheme.calculation(d0,a0,v0,Mat_state.fint,mass_mtx,...
            damp_mtx,load_s(:,1),load_s(:,2),STEP);

    [~,GT0]=apply_conditions(1,STEP,matrix,GT0);

    
    % Solve Newton Raphson
    RTOL=SOLVER.r_tolerance(BLCK);
    DTOL=SOLVER.d_tolerance(BLCK);      
    NR1=SOLVER.NR(BLCK);
    
    error_nr=zeros(SOLVER.NR_iterations(BLCK),1);
    iter=1;
    nGT0=norm(GT0);
    GT=GT0;
    
    while iter<SOLVER.NR_iterations(BLCK)
        
        iter=iter+1;
        
        % A. Delta u calculation
        [Delta_u]=Time_Scheme.solver_1(InvK,GT,a0,v0,STEP);
        
        % A.1 Check
        if isnan(du(:,iter))
            disp('Nan in displacements')
            STEP.FAIL=1;
            break;
        end
        
        A_conver=0;
        a=1;
        while A_conver==0
            
            % 1. Calculate displacements
            %--------------------------------------------------------------
            du(:,iter)=a*Delta_u;
            d0(:,1)=d0(:,1)+du(:,iter);
            
            % 2. Update
            %--------------------------------------------------------------
            % 2.1 Deformation gradient
            [Mat_state,MAT_POINT]=update_strain(d0,Mat_state,MAT_POINT,STEP);
            
            % 2.2 Stress and Internal Forces
            [~,Int_var,Mat_state,STEP]=...
                Constitutive.update(0,STEP,Int_var,Mat_state,MAT_POINT);
            
            % 2.3 Residuum
            [GT]=Time_Scheme.calculation(d0,a0,v0,Mat_state.fint,mass_mtx,...
                damp_mtx,load_s(:,1),load_s(:,2),STEP);
            [~,GT]=apply_conditions(2,STEP,matrix,GT);
            
            % 2. Convergence
            %--------------------------------------------------------------
            if isnan(GT)   %|| NORMErec(iter,1)>emax
                error('Fallo en el NR global \n');
                STEP.FAIL=1;
            else
                [CONVER,error_nr,STEP.FAIL]=LIB.convergence(GT,nGT0,error_nr,...
                    RTOL,iter,SOLVER.NR_iterations(BLCK),STEP.FAIL);
                if CONVER==1     
                    break;
                elseif iter>6 && norm(du(:,iter))/norm(d0(:,1)) < DTOL
                    CONVER=1;
                    break;
                elseif iter>6 && ...
                        norm(du(:,iter)) < 10*eps
                    CONVER=1;
                    break;
                elseif (error_nr(iter)-error_nr(iter-1))>1e-8 && iter>3
                    d0(:,1)=d0(:,1)-du(:,iter);
                    [a,CONVER]=LIB.a_factor_NR(a,error_nr,RTOL,iter);
                    if CONVER==1     
                        break;
                    end
                else
                    A_conver=1;
                end
            end
        end
        
        % B. Update
        % B.1 Position
        for j=1:GEOMETRY.nodes
            for i=1:GEOMETRY.sp
                x_a(j,i)=GEOMETRY.x_0(j,i)+d0((j-1)*GEOMETRY.df+i,1);
            end
        end
        % B.2 Shape function
        if SOLVER.REMAPPING
%         [B,near,p,gamma_,lam_LME,REMAP,wrap,EP]=LME_EP(jacobians,...
%             volume,x_a,xg,B,near,p,gamma_,lam_LME,wrap,EP,ste);
        end
        
        % B.3 Tangent matrix
        if CONVER==0 && iter>SOLVER.NR_iterations(BLCK)/2 && NR1>1
            NR1=max(1,floor(NR1/2));
        end
        if CONVER==0 && rem(iter,NR1)==0
             if STEP.ste==1
                 [stiff_mtx,Int_var,Mat_state,STEP]=...
                 Constitutive.update(4,STEP,Int_var,Mat_state,MAT_POINT);
             else
                [stiff_mtx,Int_var,Mat_state,STEP]=...
                Constitutive.update(1,STEP,Int_var,Mat_state,MAT_POINT);
            end
            [matrix]=Time_Scheme.matrix(mass_mtx,stiff_mtx,damp_mtx,STEP);
            [InvK,~]=apply_conditions(0,STEP,matrix,0);
        elseif CONVER==1
            break;
        end
    end
    Disp_field.x_a  = x_a;
    Disp_field.d    = d0;
end
