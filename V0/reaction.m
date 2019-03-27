function [OUTPUT]=reaction(residual,OUTPUT)

    global GEOMETRY

    for i=1:OUTPUT.number
        if OUTPUT.type(i,1)==2
            R=0;
            for j=1:GEOMETRY.df*GEOMETRY.nodes
                if OUTPUT.ref_list(j,i)==1
                    R=R+residual(j);
                end
            end
            OUTPUT.inst(i)=R;
        end
    end

end



