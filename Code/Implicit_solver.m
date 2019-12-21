
function Implicit_solver(BLCK,ste,ste_p,MAT_POINT,Disp_field,Int_var,...
                Mat_state,GLOBAL,stiff_mtx,MATRIX,load_s)


    %--------------------------------------------------------------------------
    % Initialize variables
    %--------------------------------------------------------------------------
    global GEOMETRY SOLVER TIME
    
    %--------------------------------------------------------------------------
    % Initial matrices
    %--------------------------------------------------------------------------
    if SOLVER.DYN(BLCK)==1
        [MATRIX]=MATRIX.matrices(Mat_state,MAT_POINT,Disp_field.d,MATRIX,BLCK);
    end
         
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % LOOP
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ste=ste+1:SOLVER.step_final

        % 1. Forces
        load_s(:,2)=load_s(:,1);
        [load_s(:,1),out_list1]=calculate_forces...
            (ste,MAT_POINT,Disp_field,Mat_state,MATRIX,BLCK);
       
        % --------------------------------------------------------
        % 2. Predictor      
        [Disp_field,Mat_state,MAT_POINT]=implicit_predictor...
            (ste,stiff_mtx,MATRIX.mass,MATRIX.damp,load_s,Disp_field,...
            MAT_POINT,Mat_state,Int_var,BLCK);
    
        % --------------------------------------------------------

        % 3. Recompute mass and damping matrices
        if SOLVER.DYN(BLCK)==1
            [MATRIX]=MATRIX.matrices(Mat_state,MAT_POINT,Disp_field.d,MATRIX,BLCK);
        end
        
        % 4. Constitutive & Stiffness_mat
        [stiff_mtx,Int_var,Mat_state]=...
                Constitutive(2,ste,Int_var,Mat_state,MAT_POINT,BLCK);
        
        
        % 5. Storage
        if rem(ste,SOLVER.SAVE_I)==0

            ste_p=ste_p+1;
            fprintf('ste_p %i \n',ste_p);
          
            [out_list2]=LIB.reaction(Mat_state.fint);
            GLOBAL.OutputList(ste_p,:)=out_list1+out_list2;
               
            GLOBAL.d(:,ste_p)   = Disp_field.d(:,1);
            GLOBAL.a(:,ste_p)   = Disp_field.a(:,1);
            GLOBAL.v(:,ste_p)   = Disp_field.v(:,1);
            
            GLOBAL.gamma(:,ste_p)   = Int_var.gamma(:,1);
            GLOBAL.epsv(:,ste_p)    = Int_var.epsv(:,1);
            GLOBAL.dgamma(:,ste_p)  = Int_var.dgamma(:,1);
            GLOBAL.H(:,ste_p)       = Int_var.H(:,1);
            GLOBAL.eta(:,ste_p)     = Int_var.eta(:,1);
            GLOBAL.Sy(:,ste_p)      = Int_var.Sy(:,1);
            GLOBAL.Sy_r(:,ste_p)    = Int_var.Sy_r(:,1);
            GLOBAL.H(:,ste_p)       = Int_var.H(:,1);
            GLOBAL.eta(:,ste_p)     = Int_var.eta(:,1);

            
            GLOBAL.xg(:,ste_p)      = LIB.reshape_S2list(MAT_POINT,'xg');

            [GLOBAL.gamma_nds(:,ste_p)]=LIB.Ep2Ep_n...
                (GLOBAL.gamma,MAT_POINT,ste_p);

            for e=1:GEOMETRY.mat_points
                [GLOBAL.Ps(e,ste_p),GLOBAL.Qs(e,ste_p)]=...
                    LIB.invar(Mat_state.Sigma(:,1),e);   %PRESSURE
            end         
            [GLOBAL.Es(:,ste_p),GLOBAL.Es_p(:,ste_p)]=...
                LIB.strains(Mat_state.F,Mat_state.Be);
            [GLOBAL.J(:,ste_p)]=LIB.S2list(MAT_POINT,'J'); %JACOBIAN
            GLOBAL.Sigma(:,ste_p)   = Mat_state.Sigma(:,1);    %STRESS
            GLOBAL.F(:,ste_p)       = Mat_state.F(:,1);  %DEFORMATION GRADIENT
            GLOBAL.Be(:,ste_p)      = Mat_state.Be(:,1); %LEFT CAUCHY GREEN 
            if SOLVER.UW==1
                GLOBAL.pw(:,ste_p) = Mat_state.pw(:,1);
                GLOBAL.Fw(:,ste_p) = Mat_state.Fw(:,1);
            elseif SOLVER.UW==2
                GLOBAL.pw(:,ste_p) = Mat_state.pw(:,1);
                GLOBAL.dpw(:,ste_p) = Mat_state.dpw(:,1);
            end

            GLOBAL.tp(ste_p,1)=TIME{BLCK}.t(ste);
            GLOBAL.ste_p=ste_p;
      
        end

        % 6. Save info
        if ((rem(ste/SOLVER.SAVE_I,SOLVER.SAVE_F)==0) ...
                || (ste==SOLVER.step_final(BLCK)) || (SOLVER.FAIL==1))
            save(SOLVER.Output(BLCK),'MAT_POINT','GLOBAL','-append')
            if SOLVER.FAIL==1
                stop
            end
            
        end
        
        % 7. Update
        Disp_field.d(:,2)=Disp_field.d(:,1);
        Disp_field.v(:,2)=Disp_field.v(:,1);
        Disp_field.a(:,2)=Disp_field.a(:,1);
        
        Mat_state.Sigma(:,2)=Mat_state.Sigma(:,1);                
        Mat_state.F(:,2)=Mat_state.F(:,1);
        Mat_state.Be(:,2)=Mat_state.Be(:,1);
        Mat_state.fint(:,2)=Mat_state.fint(:,1);
        if SOLVER.UW==1
        	Mat_state.pw(:,2)=Mat_state.pw(:,1);
            Mat_state.Fw(:,2)=Mat_state.Fw(:,1);
        elseif SOLVER.UW==2
            Mat_state.pw(:,2)=Mat_state.pw(:,1);
            Mat_state.dpw(:,2)=Mat_state.dpw(:,1);
        end
        
        Int_var.gamma(:,2)  = Int_var.gamma(:,1);
        Int_var.epsv(:,2)   = Int_var.epsv(:,1);
        Int_var.H(:,2)      = Int_var.H(:,1);
        Int_var.eta(:,2)    = Int_var.eta(:,1);
        Int_var.Sy(:,2)     = Int_var.Sy(:,1);
        Int_var.Sy_r(:,2)   = Int_var.Sy_r(:,1);
        Int_var.H(:,2)      = Int_var.H(:,1);
        Int_var.eta(:,2)    = Int_var.eta(:,1);
        Int_var.dgamma(:,2) = Int_var.dgamma(:,1);

     end
     
end

