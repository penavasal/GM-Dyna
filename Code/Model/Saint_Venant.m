function [A,T,W]=Saint_Venant(Kt,e,Ee,BLCK)
    
    % St. Venant Material

    global MATERIAL GEOMETRY SOLVER
    
    Material=GEOMETRY.material;
    MAT=MATERIAL(BLCK).MAT;
    
    G  = MAT{4,Material(e)};
    Lam= MAT{5,Material(e)};  
    %P0 = MAT{25,Material(e)};%P0;
    I=eye(3);

    A=0;
    W=0;
    
    I1=trace(Ee);
    
    T=(P0+Lam*I1)*I+2*G*Ee;
    
    if Kt==1 || Kt==2 || Kt==4
        
           A=[Lam+2*G    Lam     0   Lam ;
                Lam    Lam+2*G   0   Lam ;
                0          0     G    0  ;
                Lam      Lam     0   Lam+2*G]; 
            
    end
    
    if SOLVER.FRAC>0
        sum=0;
        for i=1:3
            for j=1:3
                sum=sum+Ee(i,j)*Ee(i,j);
            end
        end
        W=Lam/2*I1^2 + G*sum;
        
    end
    
end