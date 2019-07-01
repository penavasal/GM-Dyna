
function [A,Sc,gamma,dgamma,zetamax,etaB,H,Be]=...
    PZ(Kt,ste,e,gamma,dgamma_,zetamax,etaB,H,F,Fold,Be,P0)


    global MATERIAL
    
    MODEL=MATERIAL.MODEL;
    Mat=MATERIAL.e;
    MAT=MATERIAL.MAT;
        
    
    % Initial values
    
    Fincr=F/Fold;
    
    % Compute Trial left cauchy-Green
    BeTr = Fincr*Be*Fincr';
    
    if isnan(BeTr)
        fprintf('Error in Green-Lagrange tensor of elem e %i \n',e);
    end
    %Ee_0  = logm(Be)/2;
    Ee_tr  = logm(BeTr)/2;
    E_tot = logm(F*F')/2;
    E_0   = logm(Fold*Fold')/2;
    deps =  E_tot- E_0;

    if isnan(Ee_tr)
        error('Error in small strain tensor of elem e %i \n',e);
    elseif isreal(Ee_tr)==0
        error('Complex in small strain tensor of elem e %i \n',e);
    end
    
    Ge(13)=0;
    
    Ge(1) = MAT(29,Mat(e));%khar;
    Ge(2) = MAT(4,Mat(e));%ghar;
    Ge(3) = MAT(33,Mat(e));%alpha;
    Ge(4) = MAT(34,Mat(e));%alphag;
    Ge(5) = MAT(19,Mat(e));%Mf;
    Ge(6) = MAT(32,Mat(e));%Mg;
    Ge(7) = P0;
    Ge(8) = MAT(35,Mat(e));%beta0;
    Ge(9) = MAT(36,Mat(e));%beta1;
    Ge(10)= MAT(37,Mat(e));%H0;
    Ge(11)= MAT(38,Mat(e));%ganma;
    Ge(12)= MAT(39,Mat(e));%Hu0;
    Ge(13)= MAT(40,Mat(e));%ganmau;
    
    % Compute principal Kirchhoff tension 
    if MODEL==4.1
        [TTe,Ee,H,A,dgamma,gamma,zetamax,etaB] = ...
            PZ_backward_Euler(ste,Ge,Ee_tr,H,Kt,gamma,dgamma_,deps,zetamax,etaB);
    elseif MODEL==4.2
        [TTe,Ee,H,A,dgamma,gamma,zetamax,etaB]=...
            PZ_forward_Euler(ste,Ge,Ee_tr,H,Kt,gamma,dgamma_,deps,zetamax,etaB);
    elseif MODEL==4.3
        [TTe,Ee,H,A,dgamma,gamma,zetamax,etaB]=...
            PZ_modified_Euler(ste,Ge,Ee_tr,H,Kt,gamma,dgamma_,deps,zetamax,etaB);
    end
    
    Be = expm(2*Ee);
    Sc = TTe/det(F);

end

function [TTe,Ee,H,aep,incrlanda,defplasdes,zetamax,etaB]=...
            PZ_backward_Euler(ste,Ge,defepr,H,Kt,defplasdes,dgamma_,deps,zetamax0,etaB)
        
    alpha = Ge(3);
    alphag= Ge(4);

    % Compute volumetric and deviatoric trial strain
    [eevtrial,eestrial,theta_e,~,~,epsedev]=invar(defepr,'STRAIN');
    [dev,des]=invar(deps,'STRAIN');
    
    % Define Mg and Mf
    [Mg,Mf]=define_M(Ge,theta_e);
    etaf  = (1 + 1/alpha)*Mf;

    % Compute p, q, eta en build the Delast from the elastic law
    [De,p,q,eta]=Delast(Ge,eestrial,eevtrial);
      
    signq=sign(q);
    
    % Vectors
    [n,d]=build_vector(alpha,Mf,eta,q,signq);
    [ng,dg]=build_vector(alphag,Mg,eta,q,signq);

     discri=n(1:2)'*De*[dev;des];
        
     % H calculation    
    if discri>0
        etaB=0;  
        [Hcalcu,zetamax]=define_H(Ge,etaf,eta,p,defplasdes,zetamax0,Mg);
    elseif discri<0
        ng(1)=-abs(ng(1));
        [Hcalcu,etaB]=define_H_u(Ge,p,eta,etaB,Mg);
        zetamax=zetamax0;
    else
        Hcalcu=0;
    end
     
    % Solution initialization
    incrlanda = dgamma_;
    defplasdescalcu=defplasdes; 
    
    aa=[eevtrial;eestrial;H;incrlanda];

    
    % Solve Newton Raphson
    imax=200;
    a=1;
    TOL=1e-4;
    ite=0;
    r0=0;
    NORMErec=0; 
    while ite<imax
          
        ite=ite+1;

        % 1. Evaluate residual 
        r1=aa(1,ite)-eevtrial+aa(4,ite)*ng(1);            
        r2=aa(2,ite)-eestrial+aa(4,ite)*ng(2);
        r3=aa(3,ite)-Hcalcu;
        r4=aa(4,ite)*(aa(3,ite)+n(1:2)'*De*ng(1:2))-n(1:2)'*De*[dev;des];
        r=[r1;r2;r3;r4];                         %Vector residuos
        
        if isnan(r)   %|| NORMErec(iter,1)>emax
            error('Fallo en el PZ \n');
        else
            % 2. Check for convergence
            if ite==1
                r0=norm(r);
            else
                [CONVER,NORMErec,a,ite]=...
                    LIB.convergence(r,r0,NORMErec,TOL,ite,imax,a);
                if CONVER==1     
                    break
                end
            end
            
            % 3. Evaluate tangent matrix for NR iteration
            Dr=stimaDR(Ge,etaf,De,q,p,eta,d,dg,aa(4,ite),...
                aa(3,ite),n,ng,defplasdescalcu,dev,des,zetamax0,...
                Mg,discri,etaB,signq);
            
            if rcond(Dr)<1e-16
                disp('Small jacobian matrix of PZ return mapping');
            end
            
            % 4. Solve for displacement increment
            da = a*Dr\r;
            aa(:,ite+1)=aa(:,ite)-da;
            
            % 5. Update 
        
                % p,q,K,G and De
            [De,p,q,eta]=Delast(Ge,aa(2,ite+1),aa(1,ite+1));

                % Vectors
            [n,d]=build_vector(alpha,Mf,eta,q,signq);
            [ng,dg]=build_vector(alphag,Mg,eta,q,signq);
        
                %discri=n(1:2)'*De*[dev;des];
      
            incredefplasdes=aa(4,ite+1)*ng(2);
            defplasdescalcu=defplasdes+abs(incredefplasdes);
        
                % H calculations  
            if discri>0
                [Hcalcu,zetamax]=define_H(Ge,etaf,eta,p,defplasdescalcu,zetamax0,Mg);
            elseif discri<0
                ng(1)=-abs(ng(1));
                [Hcalcu,etaB]=define_H_u(Ge,p,eta,etaB,Mg);
            end
        
        end

    end
    
  defplasdes=defplasdescalcu; 
  eev=aa(1,ite);
  ees=aa(2,ite); 
  %es=es +des;
  %ev=ev +dev;
  H=aa(3,ite);
  incrlanda=aa(4,ite);
    
  
  % D elastoplastic, Stress and Strain
  [aep,TTe,Ee]=aep_calculation(Kt,Ge,p,q,eev,ees,epsedev,eestrial,...
     incrlanda,dgamma_,H,Mg,Mf);
  
end

function [TTe,Ee,H,aep,incrlanda,defplasdes,zetamax,etaB]=...
            PZ_forward_Euler(ste,Ge,defepr,H,Kt,defplasdes,...
            dgamma_,deps,zetamax,etaB)

    ITER=50;
        
    alpha = Ge(3);
    alphag= Ge(4);  

    % Compute volumetric and deviatoric trial strain
    [eevtrial,eestrial,theta_e,~,~,epsedev]=invar(defepr,'STRAIN');
    [dev,des]=invar(deps,'STRAIN');
    
    dev_i=dev/ITER;
    des_i=des/ITER;
    
    ees0=eestrial-des;
    eev0=eevtrial-dev;
    defplasdes_c=defplasdes;
    
    % Define Mg and Mf
    [Mg,Mf]=define_M(Ge,theta_e);
    etaf  = (1 + 1/alpha)*Mf;
    
    incrlandasum=0;
    
    for I=1:ITER
        
        eev_i=eev0+dev_i;
        ees_i=ees0+des_i;
        
        % Compute p, q, eta en build the Delast from the elastic law
        [De,p,q,eta]=Delast(Ge,ees_i,eev_i);
        signq=sign(q);
      
        % Vectors
        [n,~]=build_vector(alpha,Mf,eta,q,signq);
        discri=n(1:2)'*De*[dev_i;des_i];
        
        [ng,~]=build_vector(alphag,Mg,eta,q,signq);
        
        
        if discri>0
            % H calculation   
            [H,zetamax]=define_H(Ge,etaf,eta,p,defplasdes_c,zetamax,Mg);
        elseif discri<0
            ng(1)=-abs(ng(1));
            [H,etaB]=define_H_u(Ge,p,eta,etaB,Mg);
        end
        % Plastic multiplier
        incrlanda = discri/(H+n(1:2)'*De*ng(1:2));
        incrlandasum=incrlandasum+incrlanda;
        
        % Plastic Strain
        incredefplasvol=incrlanda*ng(1);
        incredefplasdes=incrlanda*ng(2);
        
        defplasdes_c=defplasdes_c+abs(incredefplasdes);
        
        % Elastic strain;
        eev0=eev_i-incredefplasvol;
        ees0=ees_i-incredefplasdes;
    end
    
    
    defplasdes=defplasdes_c; 
    incrlanda=incrlandasum;
    
  
  % D elastoplastic, Stress and Strain
 [aep,TTe,Ee]=aep_calculation(Kt,Ge,p,q,eev0,ees0,epsedev,eestrial,...
     incrlanda,dgamma_,H,Mg,Mf,De);
  
end


function [TTe,Ee,H,aep,incrlanda,defplasdes,zetamax,etaB]=...
            PZ_modified_Euler(ste,Ge,defepr,H,Kt,defplasdes,...
            dgamma_,deps,zetamax,etaB)

    STOL=1e-7;    
        
    hlocal=1;
    hlocalmin=1e-10;
    s=0;
    pre=0;
        
    alpha = Ge(3);
    alphag= Ge(4);  

    % Compute volumetric and deviatoric trial strain
    [eevtrial,eestrial,theta_e,~,~,epsedev]=invar(defepr,'STRAIN');
    [dev,des]=invar(deps,'STRAIN');
    
    ees0=eestrial-des;
    eev0=eevtrial-dev;
    
    defplasdes_c=defplasdes;
    
    % Define Mg and Mf
    [Mg,Mf]=define_M(Ge,theta_e);
    etaf  = (1 + 1/alpha)*Mf;
    
    while s<1 
        
        dev_i=hlocal*dev;
        des_i=hlocal*des;
        
        % Compute p, q, eta en build the Delast from the elastic law
        [De,p,q,eta]=Delast(Ge,ees0,eev0);
        signq=sign(q);

        % Vectors
        [n,~]=build_vector(alpha,Mf,eta,q,signq);
        discri=n(1:2)'*De*[dev_i;des_i];

        [ng,~]=build_vector(alphag,Mg,eta,q,signq);


        if discri>=0
            % H calculation   
            [H,zetamax]=define_H(Ge,etaf,eta,p,defplasdes_c,zetamax,Mg);
        elseif discri<0
            ng(1)=-abs(ng(1));
            [H,etaB]=define_H_u(Ge,p,eta,etaB,Mg);
        end
        % Plastic multiplier
        incrlanda = discri/(H+n(1:2)'*De*ng(1:2));

        % Plastic Strain
        incredefplasvol=incrlanda*ng(1);
        incredefplasdes=incrlanda*ng(2);
        defplasdes1=defplasdes_c+abs(incredefplasdes);

        % Elastic strain;
        eev1=eev0+dev_i-incredefplasvol;
        ees1=ees0+des_i-incredefplasdes;

        %-------------- MODIFIED------------------------------------
        % Compute p, q, eta en build the Delast from the elastic law
        [De2,p2,q2,eta2]=Delast(Ge,ees1,eev1);
        [n2,~]=build_vector(alpha,Mf,eta2,q2,signq);
        discri2=n2(1:2)'*De2*[dev_i;des_i];
        [ng2,~]=build_vector(alphag,Mg,eta2,q2,signq);

        if sign(discri)~=sign(discri)
            discri
        end
        
        if discri>=0
            % H calculation   
            [H2,zetamax]=define_H(Ge,etaf,eta2,p2,defplasdes1,zetamax,Mg);
        elseif discri<0
            ng2(1)=-abs(ng(1));
            [H2,etaB]=define_H_u(Ge,p2,eta2,etaB,Mg);
        end
        % Plastic multiplier
        incrlanda2 = discri2/(H2+n2(1:2)'*De2*ng2(1:2));

        % Plastic Strain    
        incredefplasvol2=incrlanda2*ng2(1);
        incredefplasdes2=incrlanda2*ng2(2);

        incredefelasv=dev_i-0.5*(incredefplasvol+incredefplasvol2);
        incredefplass=0.5*(incredefplasdes2+incredefplasdes);
        incredefelass=des_i-incredefplass;


        Error=0.5*(incredefplasdes-incredefplasdes2);  
        relaerror=max([1e-16,abs(Error)/(2*abs(eestrial+incredefelass))]);


        if relaerror>STOL %Decrease the step

            factor=max(0.8*sqrt(STOL/relaerror),0.1);
            hlocal=max(factor*hlocal,hlocalmin);
            pre=-1;

        else

            eev0=eev0+incredefelasv;
            ees0=ees0+incredefelass;
            defplasdes_c=defplasdes_c+abs(incredefplass);
            
            s=s+hlocal;
            factor=min(0.9*sqrt(STOL/relaerror),1.1);
            if pre==-1
                factor=min(factor,1);
                pre=0;
            end

            hlocal=factor*hlocal;
            hlocal=max(hlocal,hlocalmin);
            hlocal=min(hlocal,1-s);


        end
    end
    
    defplasdes=defplasdes_c; 
    
  
  % D elastoplastic, Stress and Strain
 [aep,TTe,Ee]=aep_calculation(Kt,Ge,p,q,eev0,ees0,epsedev,eestrial,...
     incrlanda,dgamma_,H,Mg,Mf,De);
  
end


function [n,d]=build_vector(alpha,Mf,eta,q,signq)

        if eta>0.5 
            if signq>=0
                ns  = 1;
            else
                ns  = -1;
            end
        else
            ns=0;
        end
        nv  = (1+alpha)*(Mf-eta); 
        d   = nv;
        
        det = sqrt(ns*ns+nv*nv);
        
        nt  = -1/2*q*Mf;
        
        n   = [nv ; ns; nt]/det;

end

function [De,p,q,eta]=Delast(Ge,ees,eev)

    khar  = Ge(1);
    ghar  = Ge(2);
    p0    = Ge(7);

    p=p0*exp(khar*eev+(3*ghar*khar*(ees^2)/2));
    q=-p0*3*ees*ghar*exp(khar*eev+(3*ghar*khar*(ees^2)/2));
    eta=abs(q/p);

    K=-khar*p;
    J=-khar*q;
    G=-(ghar*p+(khar*(q^2)/(3*p)));
    De=[K J;J 3*G];

end

function [H,zetamax]=define_H(Ge,etaf,eta,p,defplas,zetamax,Mg)

    alpha = Ge(3);
    beta0 = Ge(8);
    beta1 = Ge(9);
    H0    = Ge(10);
    ganma = Ge(11);
    
    zetacalcu=-p*(1-(eta/etaf))^(-1/alpha);
    zetamax=max(zetamax,zetacalcu);
    
    HDM=(zetamax/zetacalcu)^ganma;
    
    Hf=(1-(eta/etaf))^4;
    Hv=1-(eta/Mg);
    Hs=beta0*beta1*exp(-defplas*beta0);
          
    H=-H0*p*Hf*(Hv+Hs)*HDM;

end

function [Hu,etaB]=define_H_u(Ge,p,eta,etaB,Mg)

    Hu0    = Ge(12);
    ganmau = Ge(13);

    if etaB==0
        etaB=eta;
    end

    if abs(Mg/etaB)>1
        Hu=-Hu0*p*((Mg/etaB)^ganmau);   
    else
        Hu=-Hu0*p;
    end

end

function [Mg,Mf]=define_M(Ge,theta)
    if theta>pi/6
        Mf=Ge(5);
        Mg=Ge(6);
    elseif theta<-pi/6
        Mf=Ge(5)*18/(18+6*Ge(5));
        Mg=Ge(6)*18/(18+6*Ge(6));
    else
        Mf=Ge(5)*18/(18+2*Ge(5)*(1-sin(3*theta)));
        Mg=Ge(6)*18/(18+2*Ge(6)*(1-sin(3*theta)));
    end
end

function Dr=stimaDR(Ge,etaf,De,q,p,eta,d,dg,incrlanda,...
    H,n,ng,defplasdes,dev,des,chimax,Mg,discri,etaB,signq)

    % Initial parameters
    khar  = Ge(1);
    ghar  = Ge(2);
    alpha = Ge(3);
    alphag= Ge(4);
    ganma = Ge(11);
    
    beta0 = Ge(8);
    beta1 = Ge(9);
    H0    = Ge(10);
    
    Hu0    = Ge(12);
    ganmau = Ge(13);
    
    p=-p;

    % Auxiliar derivatives
    
    % De
    dDedp=[khar 0; 0 3*ghar-(khar*(q^2)/(p^2))];
    dDedq=[0 khar; khar khar*2*abs(q)/p];

    % N and Ng
    dddp=(1+alpha)*(q/(p^2)); dddq=(1+alpha)*(-1/p);
    dnvdp=(dddp/(1+d^2)^0.5)-(dddp*d^2/(1+d^2)^1.5);
    dnvdq=(dddq/(1+d^2)^0.5)-(dddq*d^2/(1+d^2)^1.5);
    dnsdp=-signq*(dddp*d/(1+d^2)^1.5);
    dnsdq=-signq*(dddq*d/(1+d^2)^1.5);

    ddgdp=(1+alphag)*(q/(p^2)); 
    ddgdq=(1+alphag)*(-1/p);
    dngvdp=(ddgdp/(1+dg^2)^0.5)-(ddgdp*dg^2/(1+dg^2)^1.5);
    dngvdq=(ddgdq/(1+dg^2)^0.5)-(ddgdq*dg^2/(1+dg^2)^1.5);
    dngsdp=-signq*(ddgdp*dg/(1+dg^2)^1.5);
    dngsdq=-signq*(ddgdq*dg/(1+dg^2)^1.5);

    % P and Q
    dpdeev=khar*p;
    dpdees=khar*q;
    dqdeev=khar*q;
    dqdees=3*ghar*p+(khar*(q^2)/p);
    
    
    if discri>0
        % H
        Hf=(1-(eta/etaf))^4;
        Hv=1-(eta/Mg);
        Hs=beta0*beta1*exp(-beta0*defplasdes);
        chi=p*(1-eta/etaf)^(-1/alpha);
        HDM=(chimax/chi)^ganma;

        dHfdp=(4*q*(1 - q/(etaf*p))^3)/(etaf*p^2);
        dHvdp=q/(Mg*p^2);
        dHsdp=-dngsdp*beta0*Hs;
        dchidp=(1 - q/(etaf*p))^(-1/alpha) - ...
            ( q*(1 - q/(etaf*p))^(-1 - 1/alpha))/(etaf*p*alpha);
        dHDMdp=ganma*(chimax/chi)^(ganma-1)*(-chimax/chi^2)*dchidp;

        dHfdq=-(4*(1 - eta/etaf)^3)/(etaf*p);
        dHvdq=-(1/(Mg*p));
        dHsdq=-dngsdq*beta0*Hs;
        dchidq=(1/(alpha*etaf))*(1-(eta/etaf))^(-1 - 1/alpha);
        dHDMdq=ganma*(chimax/chi)^(ganma-1)*(-chimax/chi^2)*dchidq;

        dHdp=H0*(Hf*(Hv+Hs)*HDM + ...
            p*dHfdp*(Hv+Hs)*HDM+...
            p*Hf*(dHvdp+dHsdp)*HDM+...
            p*Hf*(Hv+Hs)*dHDMdp);

        dHdq=H0*p*(dHfdq*(Hv+Hs)*HDM+ ...
            Hf*(dHvdq+dHsdq)*HDM+...
            Hf*(Hv+Hs)*dHDMdq);
        
        dHsdincrelanda=-beta0^2*beta1*ng(2)*exp(-beta0*defplasdes);
    
        dHdincrelanda=H0*p*Hf*dHsdincrelanda*HDM;
    else
        
        dngvdp=-sign(ng(1))*dngvdp;
        dngvdq=-sign(ng(1))*dngvdq;
        % Hu
        if abs(Mg/etaB)>1
            dHdp=Hu0*((Mg/etaB)^ganmau);   
        else
            dHdp=Hu0;
        end
        dHdincrelanda=0;
        dHdq=0;
    end
    


    dHdeev=dHdp*dpdeev+dHdq*dqdeev;
    dHdees=dHdp*dpdees+dHdq*dqdees;

    % Components of Dr
  
    Dr11=1+incrlanda*(dngvdp*dpdeev+dngvdq*dqdeev);
    Dr12=incrlanda*(dngvdp*dpdees+dngvdq*dqdees);
    Dr13=0;
    Dr14=ng(1);

    Dr21=incrlanda*(dngsdp*dpdeev+dngsdq*dqdeev);
    Dr22=1+incrlanda*(dngsdp*dpdees+dngsdq*dqdees);
    Dr23=0;
    Dr24=ng(2);

    Dr31=-dHdeev;
    Dr32=-dHdees;
    Dr33=1;
    
    if discri>0
        Dr34=H0*p*Hf*Hs*beta0*ng(2)*HDM;
    else
        Dr34=0;
    end

    Dr41=incrlanda*(dHdeev+([dnvdp dnvdq;dnsdp dnsdq]*[dpdeev;dqdeev])'*De*ng(1:2)+...
    n(1:2)'*(dpdeev*dDedp+dqdeev*dDedq)*ng(1:2)+...
    n(1:2)'*De*[dngvdp dngvdq;dngsdp dngsdq]*[dpdeev;dqdeev])-...
    ...
    (([dnvdp dnvdq;dnsdp dnsdq]*[dpdeev;dqdeev])'*De*[dev;des]+...
    n(1:2)'*(dpdeev*dDedp+dqdeev*dDedq)*[dev;des]);

    Dr42=incrlanda*(dHdees+([dnvdp dnvdq;dnsdp dnsdq]*[dpdees;dqdees])'*De*ng(1:2)+...
    n(1:2)'*(dpdees*dDedp+dqdees*dDedq)*ng(1:2)+...
    n(1:2)'*De*[dngvdp dngvdq;dngsdp dngsdq]*[dpdees;dqdees])-...
    ...
    (([dnvdp dnvdq;dnsdp dnsdq]*[dpdees;dqdees])'*De*[dev;des]+...
    n(1:2)'*(dpdees*dDedp+dqdees*dDedq)*[dev;des]);
  
    Dr43=incrlanda;
    
    Dr44=H+n(1:2)'*De*ng(1:2) + incrlanda*dHdincrelanda;

    Dr=[Dr11 Dr12 Dr13 Dr14;Dr21 Dr22 Dr23 Dr24;
        Dr31 Dr32 Dr33 Dr34;Dr41 Dr42 Dr43 Dr44];

end
    
function [q]=s_j2(s)
    %L2-norm of deviatoric stress
    q= s(1,1)^2 + s(2,2)^2 + s(3,3)^2 +...
       s(2,1)^2 + s(1,2)^2 + ...
       s(3,1)^2 + s(1,3)^2 + ...
       s(2,3)^2 + s(3,2)^2;
    q= sqrt(q);
end

function [p,q,theta,rj2,sn,devs,rj3]=invar(s,type)

    p=s(1,1)+s(2,2)+s(3,3);
    
    devs=s-p/3*eye(3);
    
    [sn]=s_j2(devs);

    rj2=sn*sn*0.5;
    

    rj3=0;
    rj3 = rj3 + 3*devs(2,1)*devs(2,1)*(devs(1,1)+devs(2,2));
    for i=1:3
        rj3 = rj3 + devs(i,i)*devs(i,i)*devs(i,i);
    end
    rj3=rj3/3;
    
    rj23=sqrt(rj2)^3;
    if 2*rj23<1.0e-18
        sint3=0;
    else
        sint3 = -3 * sqrt(3) * rj3/2/rj23;
    end
    if sint3<-1
        sint3=-1;
    elseif sint3>1
        sint3=1;
    end
    theta = 1/3*asin(sint3);
    
    if strcmp(type,'STRAIN')
        q=sqrt(2/3)*sn;
        p=-p;
    else
        q=sqrt(3*rj2);
        p=-p/3;
    end
    
    if sint3<0
        q=-q;
    end   

end

function [BB]=rota_mat_theta(p,rj2,steff,devs,rj3)

    rj2r=abs(steff/p);
    root32=sqrt(3.0/2);
    tol=1.0e-10;
    
    %cos3t=cos(3*theta);


    BB=zeros(3,4);

    % * A1 *
    BB(1,1)=-1.0/3.0;
    BB(1,2)=-1.0/3.0;
    BB(1,3)= 0.;
    BB(1,4)=-1.0/3.0;


    BB=-BB;

    % * A2 *

    if (rj2r >tol) && (steff>tol/5)
        BB(2,3)=2*devs(2,1)*root32/steff;
        BB(2,1)=devs(1,1)*root32/steff;
        BB(2,2)=devs(2,2)*root32/steff;
        BB(2,4)=devs(3,3)*root32/steff;
    else
        BB(2,:)=[-1 -1 -1 -1];
    end

    % dJ3/d
    A3(1)=devs(2,2)*devs(3,3)+rj2/3;
    A3(2)=devs(1,1)*devs(3,3)+rj2/3;
    A3(3)=-2*devs(1,2)*devs(3,3);
    A3(4)=devs(1,1)*devs(2,2)-devs(1,2)*devs(1,2)+rj2/3;
    

    
    % * A3 *

    if (rj2 >tol)
        root3=sqrt(3.0);
        root22=2*sqrt(2);
        fact0=-root3/2.0;
        comp1=root22/steff^3;
        comp2=root3*rj3/rj2^2;
        for i=1:4
            BB(3,i)=fact0*(comp1*A3(i)-comp2*BB(2,i));
        end

    end
end

function [aep,T,E_elast]=aep_calculation(Kt,Ge,P,Q,...
    epsev,epses,epsedev,n_epsedev,dgamma,dgamma_,H,Mg,Mf,De)
        
    global SOLVER
    I=eye(3);
    
    alpha = Ge(3);
    alphag= Ge(4);
    p0    = Ge(7);
    
    [De2,P2,Q2]=Delast(Ge,epses,epsev);
    %[ng,dg]=build_vector(alphag,Mg,eta,Q);
    

    dir_d =sqrt(2/3)* epsedev / n_epsedev;
    
    % Compute principal Kirchhoff tension
    T = P*I+ sqrt(2/3)*Q*dir_d;
    
     %[p1,q1]=invar(T,'STRESS');

    % Compute principal elastic strain
    E_elast = -(1/3)*epsev*I + sqrt(3/2)*epses*dir_d;
    
     [evol2,ees2]=invar(E_elast,'STRAIN');  
    
    
    % Compute algorithmic stress-strain in principal direction
    if Kt==1 || Kt==2 || Kt==4
        
        % Compute matrix De and vector ng
        [p,q,~,rj2,steff,devs,rj3]=invar(T,'STRESS');
        eta=abs(q/p);
        signq=sign(q);
        [ng,~]=build_vector(alphag,Mg,eta,q,signq);
        [n,~]=build_vector(alpha,Mf,eta,q,signq);
        [BB]=rota_mat_theta(p,rj2,steff,devs,rj3);
        
        ngxyz=BB'*ng;
        nxyz=BB'*n;
        
        [De_xyz]=assemble(De,dir_d,Ge,epsev,epses);
        if Kt==4
            aep = De_xyz;
        else
            Dpxyz=(De_xyz*ngxyz)*(nxyz'*De_xyz)/...
                (norm(nxyz)*norm(ngxyz)*H+nxyz'*De_xyz*ngxyz);
        
            aep = De_xyz - Dpxyz;
        end
        
        [t_vec1]=LIB.E2e_in(T);
        t_vec=t_vec1-[p0;p0;0;p0];
        
        [e_vec]=LIB.E2e_in(E_elast);
        t_vec2=De_xyz*e_vec;

        t_vec-t_vec2;
    else
        aep = zeros(4,4);
    end 
    
    
end


function [A]=assemble(D,n,Ge,eev,ees)


        khar  = Ge(1);
        ghar  = Ge(2);
        p0    = Ge(7);
          
        v_exp=exp(khar*eev+(3*ghar*khar*(ees^2)/2));
        
        K=-khar*p0*v_exp;
        J=khar*3*ees*ghar*v_exp*p0;
        G1=-2*ghar*p0*v_exp;
        G2=-6*khar*v_exp*p0*ees*ees*ghar*ghar;
        
        GG=3/2*(G1+G2);
        
        r23=sqrt(2/3); 
        
        m=[1 1 0 1];
        [n_vec]=LIB.E2e(n);
        
        aux=n_vec(4);
        n_vec(4)=n_vec(3);
        n_vec(3)=aux;
        
        I4=eye(4);
        I1=m'*m;
        
        M2=m'*n_vec';
        M3=n_vec*m;
        M4=n_vec*n_vec';
        
        M1=I4-1/3*I1;
              
        A2 = D(1,1)*I1 + r23*D(1,2)*M2 + r23*D(2,1)*M3 + ...
            (2/3*D(2,2)-G1)*M4 + G1*M1;
        
        %A = D(1,1)*I1 - r23*D(1,2)*M2 - r23*D(2,1)*M3 + ...
        %    2/3*D(2,2)*M1;
        
        %A = D(1,1)*I1 + r23*D(1,2)*M2 + r23*D(2,1)*M3 + ...
        %    2/3*D(2,2)*I4;

        A = K*I1 + r23*J*(M2+M3) + ...
            G2*M4 + G1*M1;
        
        p0;
        
end