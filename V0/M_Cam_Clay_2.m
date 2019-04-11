
function [A,Sc,gamma,dgamma_,Pc,Pd,Be]=...
    M_Cam_Clay_2(Kt,ste,e,gamma,dgamma_,Pcn,Pcd,F,Fold,Be,P0)

    global MODEL MAT Material
    
    % Initial values
    Mat=Material;
    
    Fincr=F/Fold;
    
    % Compute Trial left cauchy-Green
    BeTr = Fincr*Be*Fincr';
    
    if isnan(BeTr)
        fprintf('Error in Green-Lagrange tensor of elem e %i \n',e);
    end
    Ee_tr = logm(BeTr)/2;
    if isnan(Ee_tr)
        fprintf('Error in small strain tensor of elem e %i \n',e);
    end

    % Compute principal Kirchhoff tension 
    Ge(10)=0;
    
    Ge(2) = MAT(4,Mat(e));   %mu0
    Ge(3) = MAT(20,Mat(e));  %alfa
    Ge(4) = MAT(22,Mat(e));  %kappa
    Ge(5) = MAT(21,Mat(e));  %lambda
    Ge(6) = MAT(19,Mat(e));  %M
    Ge(9) = MAT(23,Mat(e));  %epsev0
    Ge(10)= MAT(30,Mat(e));  %mu
    

    if MODEL(Mat(e))==3.0
        [TTe,Ee,Pc,A,P,Q,dgamma] = tensCC(Ge,Ee_tr,Pcn,Kt,P0,dgamma_);
        Pd=Pc;
    elseif MODEL(Mat(e))==3.1
        [TTe,Ee,Pc,Pd,A,P,Q,dgamma] = visco(Ge,Ee_tr,Pcn,Pcd,Kt,P0,dgamma_,ste);
    end
    
    gamma=gamma-dgamma*2*Q/P/Ge(6)^2;
    
    Be = expm(2*Ee);
    Sc = TTe/det(F);

end


function [tenspr,epse,Pc,aep,P,Q,dgamma] = tensCC(Ge,defepr,Pcn,Kt,P0,dgamma_)
%--------------------------------------------------------
% tensVM: 
%   Compute the Kirchhoff principal tension according to CC yield criterion.
%
% Syntax:
%   [tenspr,epse,Pc,aep,flag] = tensCC(Ge,defepr,Pcn)
%
% Input:
% Ge      : Material property.
% defepr  : Elastic principal deformation.
% Pcn     : Preconsolidation pressure at step n.
%
% Output: 
% tenspr : Vector of Kirchhoff principal tension. [beta(1) beta(2) beta (3)]' 
%   epse : Vector of elastic principal strain. [epse(1) epse(2) epse (3)]'
%   Pc   : Preconsolidation pressure at step n+1.
%  aep   : Algorithmic stress-strain matrix in pr direction (a_AB).
%  flag  : Flag for plasticity.
%
% Reference:
% Borja R.I., Tamagnini C., Cam-Clay plasticity, part III: Estension of the
% infinitesimal model to include finite strain, CMAME 155, (1998), 73-95.
%
% Date:
%   Version 1.0   20.04.2018
%
% Created by: Pedro Navas
%--------------------------------------------------------
    global INIT_file

    I = eye(3);

    imax = 30;
    toll = 10e-4;

    % Set isotropic elasto-plastic parameter
    mu0  = Ge(2);  
    alfa = Ge(3);
    kappa = Ge(4);
    lambda = Ge(5);
    M = Ge(6);
    %P0 = Ge(7);
    epsev0 = Ge(9);


    % Compute volumetric trial strain
    epsevTR = defepr(1,1)+defepr(2,2)+defepr(3,3);

    % Compute deviatoric trial strain
    epsedev = defepr-1/3*I*epsevTR;

    epsesTR = sqrt(2/3)*s_j2(epsedev);


    % Compute trial stress invariants
    [Ptr,Qtr] = PQ( epsevTR, epsesTR, P0, alfa, kappa, epsev0, mu0);

    % Check for plasticity    
    Ftr = (Qtr/M)^2+Ptr*(Ptr-Pcn);

    if Ftr > toll % PLASTIC STEP    
       % Solve NR system
        % Inizialize variables
        x = zeros(3,1);
        x(1,1) = epsevTR;   % epsev
        x(2,1) = epsesTR;   % epses
        x(3,1) = 0;         % dgamma

        Pc = Pcn;
        P = Ptr;
        Q = Qtr;

        for iter = 1:imax
            % evaluate residual            
            r = [x(1,1) - epsevTR + x(3,1)*(1-(Q/M/P)^2);
                 x(2,1) - epsesTR + x(3,1)*(2*Q/P/M^2);
                 (Q/M)^2/P + P-Pc];

            if iter == 1
                r0 = r;
            end

            NORMErec(iter,1) = norm(r)/norm(r0);  
            
            [CONVER]=convergence(r,r0,NORMErec,toll,imax);

            % check for convergence
            if CONVER==1     
                break
            else
                % evaluate tangent matrix for NR iteration
                A = Atang(x(1,1),x(2,1),x(3,1),P,Q,Pc,lambda,kappa,mu0,alfa,P0,epsev0,M);
                
                if rcond(A)<1e-15
                    %disp('Small jacobian matrix of Modified Cam Clay return mapping');
                end
                % solve for displacement increment
                dx = - (A\r);
                x = x + dx;                 

                % Update 
                [P,Q]=PQ(x(1,1),x(2,1),P0,alfa,kappa,epsev0,mu0);

                OMEGA = 1/(lambda-kappa);
                Pc = Pcn*exp(-OMEGA*(epsevTR-x(1,1)));

                if iter == imax
                    if std(norm(r0)*NORMErec(iter-10:iter-1))<toll*10
                        break;
                    else
                        fprintf('\n No convergence RM \n')
                        stop;
                    end
                end
            end
        end

        %Update variable
        epsev = x(1,1);
        epses = x(2,1);
        dgamma = x(3,1);

    else % ELASTIC STEP
        %Update variable
        epsev = epsevTR;
        epses = epsesTR;

        Pc = Pcn;
        P = Ptr;
        Q = Qtr;
        
        dgamma = 0;
    end

    % Compute vector n
    n_epsedev=s_j2(epsedev);
    if n_epsedev ==0
        n = [1/sqrt(3) 0 0;
            0 1/sqrt(3) 0;
            0 0 1/sqrt(3)]; 
    else
        n = epsedev / n_epsedev;
    end

    % Compute principal Kirchhoff tension
    tenspr = P*I+ sqrt(2/3)*Q*n;

    % Compute principal elastic strain
    epse = (1/3)*epsev*I + sqrt(3/2)*epses*n;

    % Compute algorithmic stress-strain in principal direction
    if Kt==1 || Kt==2
        
        if INIT_file==0 
            dg=dgamma;
        else
            dg=dgamma_;
        end
        % Compute matrix Dep (Eq. 3.50 - Borja & Tamagnini)
        Dep = DEPtens(epsev,epses,dg,P,Q,Pc,lambda,kappa,mu0,alfa,P0,epsev0,M);
               
        if Q==0
            t1=0;
        else
            t1=2*Q/(3*epsesTR);
        end
        
        r23=sqrt(2/3); 
        
        m=[1 1 0 1];
        [n_vec]=AUX.E2e(n);
        
        I4=eye(4);
        I1=m'*m;
        
        M2=m'*n_vec';
        M3=n_vec*m;
        M4=n_vec*n_vec';
        
        M1=I4-1/3*I1;
        %M1(3,3)=0.5*M1(3,3);              % Shear term
              
        aep = Dep(1,1)*I1 + r23*Dep(1,2)*M2 + r23*Dep(2,1)*M3 + ...
            2/3*Dep(2,2)*M4 + t1*(M1-M4);
    else
        aep = zeros(4,4);
    end
    
    Kt;
end

function [tenspr,epse,Pc,Pd,aep,P,Q,dgamma] = visco(Ge,defepr,Pcn,Pcd,Kt,P0,dgamma_,ste)
    global time time_step INIT_file
    
    if ste==1
        delta_t=time_step;
    else
        delta_t=time(ste)-time(ste-1);
    end
    
    I = eye(3);

    imax = 200;
    toll = 10e-2;

    % Set isotropic elasto-plastic parameter
    mu0  = Ge(2);  
    alfa = Ge(3);
    kappa = Ge(4);
    lambda = Ge(5);
    M = Ge(6);
    %P0 = Ge(7);
    epsev0 = Ge(9);
    mu = Ge(10);
    
    tau = 24*60*60;
    
    OMEGA = 1/(lambda-kappa);
    PSI = 1/mu/OMEGA;
    XI = delta_t*mu/tau;


    % Compute volumetric trial strain
    epsevTR = defepr(1,1)+defepr(2,2)+defepr(3,3);

    % Compute deviatoric trial strain
    epsedev = defepr-1/3*I*epsevTR;

    epsesTR = sqrt(2/3)*s_j2(epsedev);


    % Compute trial stress invariants
    [Ptr,Qtr] = PQ( epsevTR, epsesTR, P0, alfa, kappa, epsev0, mu0);
    
    % Check for plasticity    
    Ftr = (Qtr/M)^2+Ptr*(Ptr-Pcn);
    
    if Ftr > toll % PLASTIC STEP    
       % Solve NR system
        % Inizialize variables
        x = zeros(4,1);
        x(1,1) = epsevTR;   % epsev
        x(2,1) = epsesTR;   % epses
        x(3,1) = 0;         % dgamma
        x(4,1) = Pcn;       % preconsolidation

        Pd = Pcd;
        P  = Ptr;
        Q  = Qtr;

        for iter = 1:imax
            
            
            % evaluate residual            
            r(:,iter) = [x(1,iter) - epsevTR + x(3,iter)*(1-(Q/M/P)^2);
                 x(2,iter) - epsesTR + x(3,iter)*(2*Q/P/M^2);
                 (Q/M)^2/P + P-x(4,iter);
                 x(3,iter)*(1-(Q/M/P)^2) + XI*(x(4,iter)/Pd)^PSI ];


            if iter>8
                r;
            end

            NORMErec(iter,1) = norm(r(:,iter))/norm(r(:,1));  
            
            [CONVER]=convergence(r(:,iter),r(:,1),NORMErec,toll,imax);

            % check for convergence
            if CONVER==1     
                break
            else
                % evaluate tangent matrix for NR iteration
                A = Atang_v(x(1,iter),x(2,iter),x(3,iter),P,Q,x(4,iter),Pd,lambda,kappa,mu,...
                    mu0,alfa,P0,epsev0,M,delta_t,tau);  
                if rcond(A)<1e-15
                    %disp('Small jacobian matrix of Modified Cam Clay return mapping');
                end
                % solve for displacement increment
                dx = - (A\r(:,iter));
                x(:,iter+1) = x(:,iter) + dx;                 

                % Update 
                [P,Q]=PQ(x(1,iter+1),x(2,iter+1),P0,alfa,kappa,epsev0,mu0);

                Pd = Pcd*exp(-OMEGA*(epsevTR-x(1,iter+1)));

                if iter == imax
                    if std(norm(r0)*NORMErec(iter-10:iter-1))<toll*10
                        break;
                    else
                        fprintf('\n No convergence RM \n')
                        stop;
                    end
                end
            end
        end

        %Update variable
        epsev  = x(1,iter);
        epses  = x(2,iter);
        dgamma = x(3,iter);
        Pc     = x(4,iter);

    else % ELASTIC STEP
        %Update variable
        epsev = epsevTR;
        epses = epsesTR;

        Pc = Pcn;
        P  = Ptr;
        Q  = Qtr;
        Pd = Pcd;
        
        dgamma = 0;
    end

    % Compute vector n
    n_epsedev=s_j2(epsedev);
    if n_epsedev ==0
        n = [1/sqrt(3) 0 0;
            0 1/sqrt(3) 0;
            0 0 1/sqrt(3)]; 
    else
        n = epsedev / n_epsedev;
    end

    % Compute principal Kirchhoff tension
    tenspr = P*I+ sqrt(2/3)*Q*n;

    % Compute principal elastic strain
    epse = (1/3)*epsev*I + sqrt(3/2)*epses*n;

    % Compute algorithmic stress-strain in principal direction
    if Kt==1 || Kt==2
        
        if INIT_file==0 
            dg=dgamma;
        else
            dg=dgamma_;
        end
        % Compute matrix Dep (Eq. 3.50 - Borja & Tamagnini)
        Dep = DEPtens(epsev,epses,dg,P,Q,Pc,lambda,kappa,mu0,alfa,P0,epsev0,M);
               
        if Q==0
            t1=0;
        else
            t1=2*Q/(3*epsesTR);
        end
        
        r23=sqrt(2/3); 
        
        m=[1 1 0 1];
        [n_vec]=AUX.E2e(n);
        
        I4=eye(4);
        I1=m'*m;
        
        M2=m'*n_vec';
        M3=n_vec*m;
        M4=n_vec*n_vec';
        
        M1=I4-1/3*I1;
        %M1(3,3)=0.5*M1(3,3);              % Shear term
              
        aep = Dep(1,1)*I1 + r23*Dep(1,2)*M2 + r23*Dep(2,1)*M3 + ...
            2/3*Dep(2,2)*M4 + t1*(M1-M4);
    else
        aep = zeros(4,4);
    end
end

function [P,Q] = PQ(epsev, epses, P0, alfa, kappa, epsev0, mu0)
    OMEGA = -(epsev-epsev0)/kappa;
    P = P0*exp(OMEGA)*(1+(3*alfa)*(epses^2)/(2*kappa));
    Q = 3*(mu0-alfa*P0*exp(OMEGA))*epses;
end

function Dep = DEPtens( epsev,epses,dgamma,P,Q,Pc,lambda,kappa,mu0,alfa,P0,epsev0,M )
    %--------------------------------------------------------
    % DEPtens: compute matrix Dep according to equatione 3.50 
    %
    % Reference:
    % Borja R.I., Tamagnini C., Cam-Clay plasticity, part III: Estension of the
    % infinitesimal model to include finite strain, CMAME 155, (1998), 73-95.
    %
    % Date: 25/10/2013
    %   Version 1.0   
    % 
    % Created by: Nicolò Spiezia
    %--------------------------------------------------------

    % Inizialize matrices
    De = zeros(2);
    % Compute parameter
    OMEGA = -(epsev-epsev0)/kappa;
    mue = mu0+(alfa/kappa)*(-P0*kappa*exp(OMEGA));

    % Compute matrix D

    De(1,1) = -P/kappa;             
    De(2,2)= 3*mue;
    De(1,2)= (3*P0*alfa*epses/kappa)*exp(OMEGA);
    De(2,1)= (3*P0*alfa*epses/kappa)*exp(OMEGA);
    
    if dgamma~=0        
        [H,b,Dp] = deal(zeros(2));  
    
        THETA = 1/(lambda-kappa);
    
        % Compute matrix H

        H(1,1) = 2;
        H(2,2) = 2/(M^2);
        H(1,2) = 0;
        H(2,1) = 0;

        % Compute matrix G

        G = H*De;

        % Compute matrix b

        b(1,1) = 1+dgamma*(G(1,1)-THETA*Pc);
        b(1,2) = dgamma*G(1,2);
        b(2,1) = dgamma*G(2,1);
        b(2,2) = 1+dgamma*G(2,2);

        % Compute parameters

        c1 = 1-dgamma*(-THETA*Pc)*(-1); 
        c2 = -dgamma*(-THETA*Pc)*0;

        d1 = De(1,1)*(2*P-Pc)+De(2,1)*(2*Q/(M^2))+THETA*Pc*(-P);
        d2 = De(1,2)*(2*P-Pc)+De(2,2)*(2*Q/(M^2));

        e = d1*(b(2,2)*(2*P-Pc)-b(1,2)*(2*Q/(M^2)))+d2*(b(1,1)*(2*Q/(M^2))-b(2,1)*(2*P-Pc));

        a1 = (d1*(b(2,2)*c1-b(1,2)*c2)+d2*(b(1,1)*c2-b(2,1)*c1)+det(b)*(-THETA*Pc)*(-P))/e;
        a2 = sqrt(2/3)*(d2*b(1,1)-d1*b(1,2))/e;

        % Compute matrix Dp

        Dp(1,1)= b(2,2)*(c1-a1*(2*P-Pc))-b(1,2)*(c2-a1*(2*Q/(M^2)));
        Dp(1,2)= b(1,2)*(-1+sqrt(3/2)*a2*(2*Q/(M^2)))-sqrt(3/2)*b(2,2)*a2*(2*P-Pc);
        Dp(2,1)= b(1,1)*(c2-a1*(2*Q/(M^2)))-b(2,1)*(c1-a1*(2*P-Pc));
        Dp(2,2)= b(1,1)*(1-sqrt(3/2)*a2*(2*Q/(M^2)))+sqrt(3/2)*b(2,1)*a2*(2*P-Pc);

        Dp = Dp/det(b);

        % Compute matrix Dep

        Dep = De*Dp;
    else
        Dep = De;
    end
end

function [A] = Atang( epsev,epses,dgamma,P,Q,Pc,lambda,kappa,mu0,alfa,P0,epsev0,M )
%--------------------------------------------------------
% ATANG Compute tangent for NR iteration for Cam Clay Return Mapping
%
% Reference:
% Borja R.I., Tamagnini C., Cam-Clay plasticity, part III: Estension of the
% infinitesimal model to include finite strain, CMAME 155, (1998), 73-95.
%
% Date: 25/10/2013
%   Version 1.0   
% 
% Created by: Nicolò Spiezia
%--------------------------------------------------------
    % Inizialize matrices
    A = zeros(3);
    [De, H] = deal(zeros(2));

    % Compute parameter
    OMEGA = -(epsev-epsev0)/kappa;
    mue = mu0+(alfa/kappa)*(-P0*kappa*exp(OMEGA));
    THETA = 1/(lambda-kappa);

    % Compute matrix D
    De(1,1) = -P/kappa;             
    De(2,2) = 3*mue;
    De(1,2) = (3*P0*alfa*epses/kappa)*exp(OMEGA);
    De(2,1) = (3*P0*alfa*epses/kappa)*exp(OMEGA);

    % Compute matrix H
    H(1,1) = 2/P*(Q/M/P)^2;
    H(2,2) = 2/P/M^2;
    H(1,2) = -2*Q/M^2/P^2;
    H(2,1) = H(1,2);

    % Compute matrix G
    G = H*De;
    
    %Derivatives
    dFdP = 1-(Q/M/P)^2;
    dFdQ = 2*Q/P/M^2;

    % Compute matrix A
    A(1,1) = 1+dgamma*G(1,1);
    A(1,2) = dgamma*G(2,1);
    A(1,3) = dFdP;

    A(2,1) = dgamma*G(1,2);
    A(2,2) = 1+dgamma*G(2,2);
    A(2,3) = dFdQ;

    A(3,1) = De(1,1)*dFdP+De(2,1)*dFdQ-THETA*Pc;
    A(3,2) = De(1,2)*dFdP+De(2,2)*dFdQ;
    A(3,3) = 0;
end

function [A] = Atang_v( epsev,epses,dgamma,P,Q,Pc,Pd,lambda,kappa,mu,...
    mu0,alfa,P0,epsev0,M,delta_t,tau)
%--------------------------------------------------------

%--------------------------------------------------------
    % Inizialize matrices
    A = zeros(4);
    [De, H] = deal(zeros(2));

    % Compute parameter
    OMEGA = -(epsev-epsev0)/kappa;
    mue = mu0+(alfa/kappa)*(-P0*kappa*exp(OMEGA));
    THETA = 1/(lambda-kappa);
    
    PSI = 1/mu/THETA;
    PSI1= PSI-1;
    XI = delta_t*mu/tau;

    % Compute matrix D
    De(1,1) = -P/kappa;             
    De(2,2) = 3*mue;
    De(1,2) = (3*P0*alfa*epses/kappa)*exp(OMEGA);
    De(2,1) = (3*P0*alfa*epses/kappa)*exp(OMEGA);

    % Compute matrix H
    H(1,1) = 2/P*(Q/M/P)^2;
    H(2,2) = 2/P/M^2;
    H(1,2) = -2*Q/M^2/P^2;
    H(2,1) = H(1,2);

    % Compute matrix G
    G = H*De;
    
    %Derivatives
    dFdP = 1-(Q/M/P)^2;
    dFdQ = 2*Q/P/M^2;

    % Compute matrix A
    A(1,1) = 1+dgamma*G(1,1);
    A(1,2) = dgamma*G(2,1);
    A(1,3) = dFdP;
    A(1,4) = 0;

    A(2,1) = dgamma*G(1,2);
    A(2,2) = 1+dgamma*G(2,2);
    A(2,3) = dFdQ;
    A(2,4) = 0;

    A(3,1) = De(1,1)*dFdP+De(2,1)*dFdQ-THETA*Pc;
    A(3,2) = De(1,2)*dFdP+De(2,2)*dFdQ;
    A(3,3) = 0;
    A(3,4) = -1;
    
    A(4,1) = dgamma*G(1,1) - THETA*XI*PSI*Pc/Pd*(Pc/Pd)^PSI1;
    A(4,2) = dgamma*G(2,1);
    A(4,3) = dFdP;
    A(4,4) = XI*PSI/Pd*(Pc/Pd)^PSI1;

end

function [q]=s_j2(s)
    %L2-norm of deviatoric stress
    q= s(1,1)^2 + s(2,2)^2 + s(3,3)^2 +...
       s(2,1)^2 + s(1,2)^2 + ...
       s(3,1)^2 + s(1,3)^2 + ...
       s(2,3)^2 + s(3,2)^2;
    q= sqrt(q);
end

function [CONVER]=convergence(r,r0,NORMErec,toll,imax)

    [iter,~]=size(NORMErec);
    CONVER=0;
    
    if norm(r) < toll 
        CONVER=1;
    elseif iter>11 && iter < imax
        long=max(round(iter/2),iter-1);
        if std(norm(r0)*NORMErec(iter-long:iter))<toll*5
            CONVER=1;
        elseif mean(NORMErec(iter-long:iter))<1e-11
            CONVER=1;
        else
            list=norm(r0)*sort(NORMErec(iter+1-long:iter));
            if std(list(1:3))<toll && std(list(long-2:long))<toll
                CONVER=1;
            elseif (std(list(1:3)) + std(list(long-2:long)))/2 < toll
                CONVER=1;
            end
        end
    elseif iter == imax
        long=max(round(iter/2),20);
        if std(norm(r0)*NORMErec(iter-long:iter))<toll*20
            CONVER=1;
        elseif mean(NORMErec(iter-long:iter))<1e-10
            CONVER=1;
        else
            list=norm(r0)*sort(NORMErec(iter+1-long:iter));
            if std(list(1:3))<toll/2 && std(list(long-2:long))<toll/2
                CONVER=1;
            elseif (std(list(1:3)) + std(list(long-2:long)))/4 < toll
                CONVER=1;
            else
                fprintf('\n No convergence RM \n')
                stop;
            end
        end
    end



end


