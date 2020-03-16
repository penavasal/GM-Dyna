function [A,T]=Neo_Hookean(Kt,e,b,J,BLCK)

    global MATERIAL GEOMETRY
    
    MODEL=MATERIAL(BLCK).MODEL;
    Material=GEOMETRY.material;
    MAT=MATERIAL(BLCK).MAT;

    if MODEL(Material(e))==1.0

        % Neo-Hookean Wriggers
        A=0;

        G  = MAT(4,Material(e));
        Lam= MAT(5,Material(e)); 
        J2=J^2;

        I=eye(3);

        T=Lam/2*(J2-1)*I+G*(b-I);
        T=T/J;

        if Kt==1 || Kt==2 || Kt==4    
               A=[Lam+2*G    J2*Lam     0       J2*Lam ;
                    J2*Lam    Lam+2*G   0       J2*Lam ;
                    0          0   (1-J2)*Lam/2+G    0  ;
                    J2*Lam   J2*Lam     0       Lam+2*G]; 
        end



    elseif MODEL(Material(e))==1.1

        % Neo-Hookean Bonet
        A=0;

        G  = MAT(4,Material(e));
        Lam= MAT(5,Material(e)); 

        I=eye(3);

        T=Lam*log(J)*I+G*(b-I);
        T=T/J;

        if Kt==1 || Kt==2 || Kt==4                   
                lam_p=Lam/J;
                mu_p=(G-Lam*log(J))/J;

                A=[lam_p+2*mu_p    lam_p       0    lam_p;
                    lam_p       lam_p+2*mu_p   0    lam_p;
                      0               0       mu_p     0;
                    lam_p          lam_p       0     lam_p+2*mu_p];
        end



    elseif MODEL(Material(e))==1.2

        % Neo-Hookean Ehlers
        A=0;

        G  = MAT( 4,Material(e));
        Lam= MAT( 5,Material(e)); 
        n0 = MAT(16,Material(e));

        if n0==0
            n0=1;
        end

        I=eye(3);

        T=Lam*n0^2*(J/n0-J/(J-1+n0))*I+G*(b-I);
        T=T/J;

        if Kt==1 || Kt==2 || Kt==4                   
            lam_p=Lam/J;
            mu_p=(G-Lam*n0^2*(J/n0-J/(J-1+n0)))/J;

            A=[lam_p+2*mu_p    lam_p       0    lam_p;
                lam_p       lam_p+2*mu_p   0    lam_p;
                  0               0       mu_p     0;
                lam_p          lam_p       0     lam_p+2*mu_p];
        end

    end

end