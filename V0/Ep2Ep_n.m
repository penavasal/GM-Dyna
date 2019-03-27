
function [Gamma_nds]=Ep2Ep_n(Gamma_tot,Shape_function,ste_p)
    
    global GEOMETRY
    Gamma_nds=zeros(GEOMETRY.nodes,1);
    
    for i=1:GEOMETRY.elements
        sh=Shape_function.p{i};
        nds=Shape_function.near{i};
        n=length(nds);
        for j=1:n
            Gamma_nds(nds(j),1)=Gamma_nds(nds(j),1)+sh(j)*Gamma_tot(i,ste_p);
        end
    end

end