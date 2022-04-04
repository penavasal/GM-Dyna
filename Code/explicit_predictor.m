% Hay que arreglar la parte de contacto

function [Disp_field,Mat_state,MAT_POINT]=explicit_predictor...
            (STEP,Disp_field,MAT_POINT,Mat_state)

    global GEOMETRY TIME SOLVER
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    nodes=GEOMETRY.nodes;
    
    time_step=STEP.dt;
    BLCK=STEP.BLCK;
    
    gamma=TIME{BLCK}.gamma;
    beta=TIME{BLCK}.beta;
    af=TIME{BLCK}.af;
    %am=TIME{BLCK}.am;
    
    
    x_a = Disp_field.x_a;
    d0  = Disp_field.d;
    a0  = Disp_field.a;
    v0  = Disp_field.v;
    
    %% 0. Boundary conditions
    [boundary,i_disp,velo]=calculate_boundaries(STEP);

    %% 1. Predictor
    C1=1;
    C2=1;
    while C1 || C2
        if af~=1
            for i=1:nodes*df
                if boundary(i)==0
                    if SOLVER.UW==2 && mod(i,3)==0
                        d0(i,1)=d0(i,2)+(1-gamma)*time_step*v0(i,2);
                    elseif SOLVER.UW==3 && mod(i,5)==0
                        d0(i,1)=d0(i,2)+(1-gamma)*time_step*v0(i,2);
                    else
                        d0(i,1)=d0(i,2)+time_step*v0(i,2)+...
                            (0.5-beta)*time_step^2*a0(i,2);
                        v0(i,1)=v0(i,2)+(1-gamma)*time_step*a0(i,2);
                    end
                elseif boundary(i)==1
                    d0(i,1)=d0(i,2)+i_disp(i,1);
                    v0(i,1)=i_disp(i,1)/time_step;
                else
                    v0(i,1)=velo(i,1);
                    d0(i,1)=d0(i,2)+velo(i,1)*time_step;
                end
            end
        end
        
        if af~=0
            daf(:,1)=d0(:,1)*(1-af)+d0(:,2)*(af);
            daf(:,2)=d0(:,3);
        else
            daf=d0;
        end
        
%         if SOLVER.UW==2
%             Mat_state=PW.predictor(Mat_state,MAT_POINT,d0);
%         end

        for j=1:nodes
            for i=1:sp
                x_a(j,i)=GEOMETRY.x_0(j,i)+d0((j-1)*df+i,1);
            end
        end

%         if MASTER
%             [a1,C1]=contact(...
%                 nCC,MASTER,M_LIST,Mat_nds,x_a1,time_step,lmass,vs,as);
%         else
%             C1=0;
%         end
%         
%         [a1,C2]=PS_surface(nPS,PS_LIST,x_a1,lmass,vs,as,x_ps,K,h_min);
          C1=0;
          C2=0;

    end
    
    %% 2. REMAPPING

         
    [MAT_POINT]=update_mp(d0,MAT_POINT);
            
    [Mat_state,MAT_POINT]=update_strain(daf,Mat_state,MAT_POINT,STEP);

    REMAP=SOLVER.REMAPPING;  %Flag
    if REMAP==1
        MAT_POINT=SH.remap(MAT_POINT,Disp_field);
        %[B,near,p,gamma_,lam_LME,REMAP,wrap,EP]=LME_EP(jacobians,...
        %             volume,x_a,xg,B,near,p,gamma_,lam_LME,wrap,EP,ste);
    end

    %% 3. Assemble vectors
    Disp_field.x_a  = x_a;
    Disp_field.d    = d0;
    Disp_field.v    = v0;
    
end
