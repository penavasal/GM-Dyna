function [A,T,Be]=Saint_Venant(Kt,e,F,BLCK)
    
    % St. Venant Material

    global MATERIAL GEOMETRY
    
    Material=GEOMETRY.material;
    MAT=MATERIAL(BLCK).MAT;
    
    G  = MAT(4,Material(e));
    Lam= MAT(5,Material(e));   
    I=eye(3);

    A=0;
    
    Be=F*F';
    Ee=logm(Be)/2;
    %Ee(2,1)=Ee(2,1)/2;
    %Ee(1,2)=Ee(1,2)/2;

    T=Lam*trace(Ee)*I+2*G*Ee;
    
    if Kt==1 || Kt==2 || Kt==4
        
           A=[Lam+2*G    Lam     0   Lam ;
                Lam    Lam+2*G   0   Lam ;
                0          0     G    0  ;
                Lam      Lam     0   Lam+2*G]; 
            
    end
    
%      E_vec(1,1)=Ee(1,1);
%      E_vec(2,1)=Ee(2,2);
%      E_vec(4,1)=Ee(3,3);
%      E_vec(3,1)=Ee(2,1);
%      T2=A*E_vec;
%      e;



end