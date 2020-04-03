function [A,T,W]=Neo_Hookean(Kt,e,b,J,P0,BLCK)

    global MATERIAL GEOMETRY SOLVER
    
    MODEL=MATERIAL(BLCK).MODEL;
    Material=GEOMETRY.material;
    MAT=MATERIAL(BLCK).MAT;

    A=0;
    W=0;
    
    I=eye(3);
    
    G  = MAT{4,Material(e)};
    Lam= MAT{5,Material(e)}; 
    
    if MODEL(Material(e))==1.0

        % Neo-Hookean Wriggers
        J2=J^2;

        T=Lam/2*(J2-1)*I+G*(b-I);
        T=T+P0(1)*I;
        T=T/J;

        if Kt==1 || Kt==2 || Kt==4    
               A=[Lam+2*G    J2*Lam     0       J2*Lam ;
                    J2*Lam    Lam+2*G   0       J2*Lam ;
                    0          0   (1-J2)*Lam/2+G    0  ;
                    J2*Lam   J2*Lam     0       Lam+2*G]; 
        end
        
        if SOLVER.FRAC>0
            
            I1=trace(b);
            K=lam+2*G/3;
            W=K/4*(J^2 - 1) + G/2*(I1-3)-K/2*log(J);

        end

    elseif MODEL(Material(e))==1.1

        % Neo-Hookean Bonet

        T=Lam*log(J)*I+G*(b-I);
        T=T+P0(1)*I;
        T=T/J;

        if Kt==1 || Kt==2 || Kt==4                   
                lam_p=Lam/J;
                mu_p=(G-Lam*log(J))/J;

                A=[lam_p+2*mu_p    lam_p       0    lam_p;
                    lam_p       lam_p+2*mu_p   0    lam_p;
                      0               0       mu_p     0;
                    lam_p          lam_p       0     lam_p+2*mu_p];
        end
        
        if SOLVER.FRAC>0
            
            I1=trace(b);
            W=Lam/2*log(J)^2 + G/2*(I1-3)-G*log(J);

        end



    elseif MODEL(Material(e))==1.2

        % Neo-Hookean Ehlers
        n0 = MAT(16,Material(e));

        if n0==0
            n0=1;
        end

        T=Lam*n0^2*(J/n0-J/(J-1+n0))*I+G*(b-I);
        T=T+P0(1)*I;
        T=T/J;

        if Kt==1 || Kt==2 || Kt==4                   
            lam_p=Lam/J;
            mu_p=(G-Lam*n0^2*(J/n0-J/(J-1+n0)))/J;

            A=[lam_p+2*mu_p    lam_p       0    lam_p;
                lam_p       lam_p+2*mu_p   0    lam_p;
                  0               0       mu_p     0;
                lam_p          lam_p       0     lam_p+2*mu_p];
        end
        
        if SOLVER.FRAC>0
            
            I1=trace(b);
            W=G/2*(I1-3)-G*log(J);
            W=W+lam*(1-n0)^2*((J-1)/(1-n0)-ln((J-n0)/(1-n0)));

        end

    end

end