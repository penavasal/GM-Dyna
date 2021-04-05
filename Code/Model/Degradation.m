
function [A,Sc,gamma,dgamma,sy,Ee,E_ini]=...
    Degradation(Kt,e,gamma,dgamma,Ee,Edev,E_ini,P0,BLCK,dt,MAT_POINT)

    global MATERIAL GEOMETRY
    
    
    Mat=GEOMETRY.material;
    MAT=MATERIAL(BLCK).MAT;
    
    % Initial values
    I=eye(3);
    A=zeros(4,4);
    
    K=MAT{5,Mat(e)}+2*MAT{4,Mat(e)}/3;  % Bulk modulus
    G = MAT{4,Mat(e)};                  % Shear modulus
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                     RADIAL RETURN MAPPING                       %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    tol=1e-11;                   % Error tolerance		
    loop=20;			        % Positive loop

    %------------------
    %Elastic predictor
    %------------------
    
    % Deviatoric and volumetric strain
    [e_dev,e_vol]=split(Ee);
    p=-(K*e_vol+P0(1));
    s=2*G*e_dev;
    
    %Sig2=2*MAT(4,Mat(e))*Ee+MAT(5,Mat(e))*trace(Ee)*I;
    %[s2,p2]=split(Sig2);
    
    
    [snorm]=s_j2(s);             %L2-norm of deviatoric stress
    q = sqrt(3/2)*snorm;
    %[snorm2]=s_j2(s2);
    %sy=-sy;
    
    mat95(1) = -MAT{62,Mat(e)};
    mat95(2) = MAT{64,Mat(e)};
    mat95(3) = -MAT{7,Mat(e)};
    [sy,H,E_ini]=FRAC.eigendegradation(e,Edev,E_ini,MAT_POINT,mat95);
    
    if (q - sy)/sy < tol
        dEp=zeros(3);
        Sc=-p*I+s;
        dg=0;
    else
    %------------------
    % Plastic corrector
    %------------------
        if E_ini(2)==0
            E_ini(2)=E_ini(1);
        end
    
        N_1=1/MAT{14,Mat(e)};
        mu=MAT{13,Mat(e)};
        
        g=gamma;            % gamma = gamma_0
        
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

                f(iter)=(q - 3*G*dgm(iter))/sy - (1+(dgm(iter)/mu/dt)^N_1);  % fi=0
                df=-H*(q - 3*G*dgm(iter))/sy^2 - 3*G/sy - N_1/mu/dt*(dgm(iter)/mu/dt)^(N_1-1);       % dfi/dgamma

                iter=iter+1;
                ddg(iter)=t*f(iter-1)/df;              %increment of dgamma
                
                dgm(iter)=dgm(iter-1)+ddg(iter);

                gm(iter)=gm(iter-1)+ddg(iter);


                count=0;
                while (dgm(iter)<0.0 && count<loop)

                    dgm(iter)=dgm(iter-1)/2;
                    gm(iter)=gm(1)+dgm(iter);
                    
                    f(iter-1)=(q - 3*G*dgm(iter))/sy - (1+(dgm(iter)/mu/dt)^N_1);  % fi=0
                    df=-H*(q - 3*G*dgm(iter))/sy^2 - 3*G/sy - N_1/mu/dt*(dgm(iter)/mu/dt)^(N_1-1);       % dfi/dgamma

                    ddg(iter)=t*f(iter-1)/df;             %increment of dgamma
                    dgm(iter)=dgm(iter-1)+ddg(iter);
                    gm(iter)=gm(iter-1)+ddg(iter);

                end

                if abs(f(iter-1))>=error
                    ddg(iter)=ddg(iter)/10;
                    dgm(iter)=dgm(iter-1)+ddg(iter);
                    gm(iter)=gm(iter-1)+ddg(iter);
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
        dg=dgm(iter);
        
        d1=1-3*G*dg/q;
        dEp=dg*s/snorm;
    
        Sc=-p*I+d1*s;
    end
    
    %sy=-sy;
%     [ss,p]=split(Sc);
%     [q]=s_j2(ss);
%     Q= sqrt(3/2)*q
%     P=p/3
    %-----------------------------------------------
    % Update tensors and plastic variables
    %-----------------------------------------------  
    if isnan(dEp) | (abs(rcond(dEp)))<1e-8 | isnan(rcond(dEp))
        e;
    end

    %Update strain
    Ee = Ee - dEp;
    
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
        [s_vec]=LIB.E2e(s);
        
        I1=m'*m;
        M1=I4-1/3*I1;
        M1(3,3)=0.5*M1(3,3);              % Shear term
        M2=m'*s_vec';
        M3=s_vec*m;
        M4=s_vec*s_vec';
        
        A=(D3*I1+D2*M1-D5*M2-D6*M3-D4*M4)*D7;
        
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
