
function [A,Sc,gamma,dgamma,sy,Be]=...
    Drucker_prager(Kt,e,gamma,dgamma,sy,F,Be,Fold)

    global MATERIAL
    
    MODEL=MATERIAL.MODEL;
    Mat=MATERIAL.e;
    MAT=MATERIAL.MAT;
    
    % Initial values
    I=eye(3);
    A=zeros(4,4);
    
    K=MAT(5,Mat(e))+2*MAT(4,Mat(e))/3;  % Bulk modulus
    G = MAT(4,Mat(e));                  % Shear modulus
    
    % Predictor tensors
    Fincr=F/Fold;
    
    % Compute Trial left cauchy-Green
    BeTr = Fincr*Be*Fincr';   

    if isnan(BeTr)
        fprintf('Error in Green-Lagrange tensor of elem e %i \n',e);
    end
    Ee = logm(BeTr)/2;
    if isnan(Ee)
        fprintf('Error in small strain tensor of elem e %i \n',e);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                     RADIAL RETURN MAPPING                       %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    tol=1e-8;                   % Error tolerance		
    loop=20;			        % Positive loop
    
    cte=sqrt(2/3);
    
    if MODEL==2.0
        alfn=0;
        alfm=0;
        beta=cte;
    elseif MODEL==2.1
        alfn=2*cte*sin(MAT(11,Mat(e)))/(3-sin(MAT(11,Mat(e))));  % Alpha_n
        alfm=2*cte*sin(MAT(12,Mat(e)))/(3-sin(MAT(12,Mat(e))));  % Alpha_m
        beta=cte*6*cos(MAT(11,Mat(e)))/(3-sin(MAT(11,Mat(e))));  % Beta*sqrt(2/3)
    elseif MODEL==2.2
        alfn=2*cte*sin(MAT(11,Mat(e)))/(3+sin(MAT(11,Mat(e))));  % Alpha_n
        alfm=2*cte*sin(MAT(12,Mat(e)))/(3+sin(MAT(12,Mat(e))));  % Alpha_m
        beta=cte*6*cos(MAT(11,Mat(e)))/(3+sin(MAT(11,Mat(e))));  % Beta*sqrt(2/3)
    elseif MODEL==2.3
        alfn=cte*tan(MAT(11,Mat(e)))/sqrt(3+4*(tan(MAT(11,Mat(e))))^2);  % Alpha_n
        alfm=cte*tan(MAT(12,Mat(e)))/sqrt(3+4*(tan(MAT(12,Mat(e))))^2);  % Alpha_m
        beta=cte*3/sqrt(3+4*(tan(MAT(11,Mat(e))))^2);  % Beta*sqrt(2/3)
    end
    
    %------------------
    %Elastic predictor
    %------------------
    
    % Deviatoric and volumetric strain
    [e_dev,e_vol]=split(Ee);
    p=K*e_vol;
    s=2*G*e_dev;
    
    %Sig2=2*MAT(4,Mat(e))*Ee+MAT(5,Mat(e))*trace(Ee)*I;
    %[s2,p2]=split(Sig2);
    
    
    [snorm]=s_j2(s);             %L2-norm of deviatoric stress
    %[snorm2]=s_j2(s2);
    
    if snorm + 3*alfn*p - beta*sy < tol
        dEp=zeros(3);
        Sc=p*I+s;
        %Sc2=p2*I+s2;
        dg=0;
        dg1=0;
    else
    %------------------
    % Plastic corrector
    %------------------
        N_1=1/MAT(9,Mat(e));
        ads = sqrt(3*alfm*alfm+1);
        dg=0;               % initial delta_gamma = 0
        g=gamma;            % gamma = gamma_0
        [H]=der_Sy(MAT(7,Mat(e)),MAT(10,Mat(e)),N_1,g);
        
        %%%%
        p_lim = (4.5*K*alfm*snorm/G + beta/alfn*(snorm*ads*H/2/G + sy))/3;
        
        if p<p_lim
            dg1=0;
            error=1e32;
            iter=1;
            gm(iter)=g;
            dgm(iter)=0;
            ddg(iter)=1e32;
            %-------------CLASSIC----ITERATOR------------------------------
            while error>tol && abs(ddg(iter))>2*eps  &&  iter<300
                if(iter<300) 
                    t=-1;
                    
                    [s_y]=sig_y(MAT(7,Mat(e)),MAT(10,Mat(e)),N_1,gm(iter));  % Sigma_y (Hardening law)

                    if s_y<0
                        s_y=0;
                        H=0;
                        %t=1;
                    else
                        [H]=der_Sy(MAT(7,Mat(e)),MAT(10,Mat(e)),N_1,gm(iter));   % H=dS/dE (Hardening law)
                    end

                    f(iter)=snorm - 2*G*dgm(iter) + 3*alfn*(p-3*K*alfm*dgm(iter)) - beta*s_y;  % fi=0
                    df=-9*K*alfn*alfm-2*G-H*beta*ads;                     % dfi/dgamma

                    iter=iter+1;
                    ddg(iter)=t*f(iter-1)/df;                     %increment of dgamma
                    dgm(iter)=dgm(iter-1)+ddg(iter);

                    gm(iter)=gm(iter-1)+ddg(iter)*ads;


                    count=0;
                    while (dgm(iter)<0.0 && count<loop)

                        dgm(iter)=dgm(iter-1)/2;
                        gm(iter)=gm(1)+dgm(iter)*ads;

                        if s_y<0
                            s_y=0;
                            H=0;
                            %t=1;
                        else
                            [H]=der_Sy(MAT(7,Mat(e)),MAT(10,Mat(e)),N_1,gm(iter));   % H=dS/dE (Hardening law)
                        end

                        f(iter-1)=snorm - 2*G*dgm(iter) + 3*alfn*(p-3*K*alfm*dg) - beta*s_y;  % fi=0
                        df=-9*K*alfn*alfm-2*G-H*beta*ads;                     % dfi/dgamma

                        ddg(iter)=t*f(iter-1)/df;                     %increment of dgamma
                        dgm(iter)=dgm(iter-1)+ddg(iter);
                        gm(iter)=gm(iter-1)+ddg(iter)*ads;

                    end
                    
                    if abs(f(iter-1))>=error
                        ddg(iter)=ddg(iter)/10;
                        dgm(iter)=dgm(iter-1)+ddg(iter);
                        gm(iter)=gm(iter-1)+ddg(iter)*ads;
                    else
                        error=abs(f(iter-1));
                    end
                else
                    disp('Problema de convergencia')
                    %return
                end
            end
            %--------------END----CLASSIC----ITERATOR-----------------------
            gamma=gm(iter);
            dgamma=dgm(iter);
            d1=1-2*G*dgamma(e,1)/snorm;
            dEp=dgamma(e,1)*alfm*I+dgamma(e,1)*s/snorm;
        else
            t=-1;
            error=1e32;
            iter=1;
            dg1 = snorm/2/G;
            dg2(iter) = 0;
            coef = sqrt(dg1*dg1+3*alfm*alfm*(dg1+dg2(iter))^2);
            ddg(iter) = 0;
            CON=0;
            
            if alfm~=0
         
                [s_y]=sig_y(MAT(7,Mat(e)),MAT(10,Mat(e)),N_1,g);  % Sigma_y (Hardening law)
                if s_y<0
                    s_y=0;
                    H=0;
                    %t=1;
                else
                    [H]=der_Sy(MAT(7,Mat(e)),MAT(10,Mat(e)),N_1,g);   % H=dS/dE (Hardening law)
                end
                %---------------APEX---ITERATOR---------------------------------
                while error>tol  &&  CON==0
                    if(iter<150)
                        
                        coef = sqrt(dg1*dg1+3*alfm*alfm*(dg1+dg2(iter))^2);

                        f(iter)=beta*H*coef/3/alfn + 3*K*alfm*dg2(iter) - ...
                            p+3*K*alfm*dg1+beta*s_y/3/alfn;                       % fi=0
                        df(iter)=3*K*alfm + 3*alfm*alfm*(dg1+dg2(iter))*H*beta/(3*coef*alfn); % dfi/dgamma

                        
                        iter=iter+1;
                        ddg(iter)=t*f(iter-1)/df(iter-1);                     %increment of dgamma
                        dg2(iter)=dg2(iter-1)+ddg(iter);

                        count=0;
                        
                        % Avoid oscillations
                        if iter>3 && ddg(iter)==ddg(iter-2)
                           if ddg(iter)==ddg(iter-1)
                               f(iter-1)=0;
                           else
                               ddg(iter)=ddg(iter)/2;
                               dg2(iter)=dg2(iter-1)+ddg(iter);
                               if dg2(iter)==dg2(iter-1)
                                   f(iter-1)=0;
                               end
                           end
                        end
                        
                        % Avoid negative values
                        while (dg2(iter)<0.0 && count<loop) && iter>2

                            dg2(iter)=dg2(iter)/2;

                            coef = sqrt(dg1*dg1+3*alfm*alfm*(dg1+dg2(iter))^2);

                            f(iter-1)=beta*H*coef/3/alfn + 3*K*alfm*dg2(iter) - ...
                                p+3*K*alfm*dg1+beta*s_y/3/alfn;                       % fi=0
                            df(iter-1)=3*K*alfm + 3*alfm*alfm*(dg1+dg2(iter))*H*beta/(3*coef*alfn); % dfi/dgamma

                            ddg(iter)=t*f(iter-1)/df(iter-1);                     %increment of dgamma
                            dg2(iter)=dg2(iter)+ddg(iter);
                        end
                        
                        % If negative at the beginning, breackage, no
                        % plastic
                        if (dg2(iter)<0) && (iter==2)
                            CON=1;
                        end

                        error=abs(f(iter-1));
                    else
                        disp('Problema de convergencia')
                        %return
                    end

                end
            
            end
           %--------------END----APEX----ITERATOR--------------------------
           if CON==1
               dg1=0;
               dg=0;
           else
               dg=dg1+dg2(iter);
           end
           dgamma=sqrt(dg1*dg1+3*alfm*alfm*dg*dg);
           gamma=gamma+dgamma;
           s_y=max(0,sig_y(MAT(7,Mat(e)),MAT(10,Mat(e)),N_1,gamma));
           d1=0;
           dEp=dg*alfm*I+dg1*s/snorm;
        end
        sy=s_y;
        Sc=(p-3*alfm*K*dg)*I+d1*s;
    end
    
    %-----------------------------------------------
    % Update tensors and plastic variables
    %-----------------------------------------------  
    if isnan(dEp) | (abs(rcond(dEp)))<1e-8 | isnan(rcond(dEp))
        e;
    end

    %Update strain
    Ee = Ee - dEp;
    Be = expm(2*Ee); 
    
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                    CONSISTENT TANGENT MATRIX                    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if Kt==1 || Kt==2 || Kt==4
        
        if Kt==4
            dg=0;
        end
        
        % Coefficients
        if dg==0            %Elastic
            D2=2*G;
            D3=K;
            D4=0;
            D5=0;
            D6=0;
            D7=1.0;
            
        elseif dg1==0       %Classic
            D1=1 - 2*G*dg/snorm;
            D2=2*G*D1;
            C2=1/(9*alfm*alfn*K + 2*G + ads*beta*H);
            D3=K-9*K^2*alfn*alfm*C2;
            D4=4*G^2*(C2-dg/snorm)/snorm/snorm;
            D5=6*K*G*alfm*C2/snorm;
            D6=6*K*G*alfn*C2/snorm;
            D7=1.0;
        else
%             D2=0;
%             D3=K;
%             D4=0;
%             D5=-K/(2*G*alfm*dg);
%             D6=0;
%             D7=alfm*beta*H*dg / (3*alfn*K*coef + alfm*beta*H*dg);
            D2=0;
            D3=K*alfm*dg;
            D4=0;
            D5=-K/(2*G);
            D6=0;
            D7=beta*H / (3*alfn*K*coef + alfm*beta*H*dg);
        end
        
        % Matrixes
        I4=eye(4);
        m=[1 1 0 1];
        [s_vec]=E2e(s);
        
        I1=m'*m;
        M1=I4-1/3*I1;
        M1(3,3)=0.5*M1(3,3);              % Shear term
        M2=m'*s_vec';
        M3=s_vec*m;
        M4=s_vec*s_vec';
        
        A=(D3*I1+D2*M1-D5*M2-D6*M3-D4*M4)*D7;

        if isnan(A)
            e;
        end
        
    end
    
end

function [A_dev,a_vol]=split(A)

    I=eye(3);
    a_vol=(A(1,1)+A(2,2)+A(3,3));
    A_dev=A-a_vol*I/3;

end

function [q]=s_j2(s)
    %L2-norm of deviatoric stress
    q= s(1,1)^2 + s(2,2)^2 + s(3,3)^2 +...
       s(2,1)^2 + s(1,2)^2 + ...
       s(3,1)^2 + s(1,3)^2 + ...
       s(2,3)^2 + s(3,2)^2;
    q= sqrt(q);
end

function [s_y]=sig_y(Sy0,E0,N_1,g)
    s_y=Sy0*(1+g/E0)^(N_1);
end
function [H]=der_Sy(Sy0,E0,N_1,g)
    H=Sy0*N_1/E0*(1+g/E0)^(N_1-1);
end

function [e]=E2e(E)
   
    e=zeros(4,1);

    %Build vector
    e(1)=E(1,1);
    e(2)=E(2,2);
    e(3)=E(1,2);
    e(4)=E(3,3);
end

        
       