% Hay que arreglar la parte de contacto

function [Disp_field,Mat_state,MAT_POINT,FAIL]=explicit_predictor...
            (ste,Disp_field,MAT_POINT,Mat_state,FAIL)

    global GEOMETRY TI_scheme TIME
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    nodes=GEOMETRY.nodes;
    
    if ste==1
        time_step=TIME.t(ste+1)-TIME.t(ste);
    else
        time_step=TIME.t(ste)-TIME.t(ste-1);
    end  
    
    gamma=TI_scheme.gamma;
    x_a = Disp_field.x_a;
    d0  = Disp_field.d;
    a0  = Disp_field.a;
    v0  = Disp_field.v;
    
    %% 0. Boundary conditions
    [boundary,i_disp,velo]=calculate_boundaries(ste);

    %% 1. Predictor
    C1=1;
    C2=1;
    while C1 || C2
        for i=1:nodes*df
            if boundary(i)==0
                d0(i,1)=d0(i,2)+time_step*v0(i,2)+0.5*time_step^2*a0(i,2);
                v0(i,1)=v0(i,2)+(1-gamma)*time_step*a0(i,2);
            elseif boundary(i)==1
                d0(i,1)=d0(i,2)+i_disp(i,1);
                v0(i,1)=i_disp(i,1)/time_step;
            else
                v0(i,1)=velo(i,1);
                d0(i,1)=d0(i,2)+velo(i,1)*time_step;
            end
        end

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
%         REMAP=1;  %Flag
%         iter=1;
%         error_tol=min(h_ini.*sqrt(jacobians))*1e-6;
%         wrap=zeros(GEOMETRY.mat_points,1);
%         while REMAP==1           
            [MAT_POINT]=update_mp(d0,MAT_POINT);
            
            [Mat_state,MAT_POINT]=update_F(d0,Mat_state,MAT_POINT);

%             if iter==1
%                 xg1=xg2;
%             else
%                 error=norm(abs(xg1-xg2));
%                 if error<error_tol
%                     REMAP=0;
%                 else
%                     xg1=xg2;
%                     if iter>=3
%                         iter;
%                     end
%                 end
%             end
% 
%             if REMAP==1
%                 [B,near,p,gamma_,lam_LME,REMAP,wrap,EP]=LME_EP(jacobians,...
%                     volume,x_a,xg,B,near,p,gamma_,lam_LME,wrap,EP,ste);
%             end
%             iter=iter+1;
%         end

    %% 3. Assemble vectors
    Disp_field.x_a  = x_a;
    Disp_field.d    = d0;
    Disp_field.v    = v0;
    
end
