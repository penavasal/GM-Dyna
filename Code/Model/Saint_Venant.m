function [A,T]=Saint_Venant(Kt,e,Ee,BLCK)
    
    % St. Venant Material

    global MATERIAL GEOMETRY
    
    Material=GEOMETRY.material;
    MAT=MATERIAL(BLCK).MAT;
    
    G  = MAT(4,Material(e));
    Lam= MAT(5,Material(e));  
    P0 = MAT(25,Material(e));%P0;
    I=eye(3);

    A=0;
    
    T=(P0+Lam*trace(Ee))*I+2*G*Ee;
    
    if Kt==1 || Kt==2 || Kt==4
        
           A=[Lam+2*G    Lam     0   Lam ;
                Lam    Lam+2*G   0   Lam ;
                0          0     G    0  ;
                Lam      Lam     0   Lam+2*G]; 
            
    end
    
end