function [es,es_p]=strains(def_G,b_e)

    global GEOMETRY

    [es,es_p]=deal(zeros(GEOMETRY.elements*4,1));
    [f_v,be]=deal(zeros(5,1));
    
    for e=1:GEOMETRY.elements
        for i=1:5
            f_v(i,1)=def_G((e-1)*5 + i,1);
            be(i,1)=b_e((e-1)*5 + i,1);
        end           
        [F]=v2m(f_v,GEOMETRY.sp);
        [Be]=v2m(be,GEOMETRY.sp);

        Btot = F*F';
        Etot = logm(Btot)/2;
        Ee   = logm(Be)/2;
        Ep   = Etot-Ee;
        
        [ee]=E2e(Ee);
        [ep]=E2e(Ep);
        for i=1:4
            es((e-1)*4+i,1)=ee(i,1);
            es_p((e-1)*4+i,1)=ep(i,1);
        end
    end

end

function [e]=E2e(E)
   
    e=zeros(4,1);

    %Build vector
    e(1)=E(1,1);
    e(2)=E(2,2);
    e(3)=E(3,3);
    e(4)=E(1,2);%+E(2,1);
end