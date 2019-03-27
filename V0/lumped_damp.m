function [C]=lumped_damp(Shape_function,Mat_state)

    global GEOMETRY SOLVER
    
    sp=GEOMETRY.sp;

    C=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
    if SOLVER.UW
        %% Lumped Damp **********************
        for i=1:GEOMETRY.elements
            volume=GEOMETRY.Area(i)*Mat_state.J(i);
            nd = Shape_function.near{i};
            m  = length(nd);
            sh = Shape_function.p{i};
            for t1=1:m
                for k=1:sp
                    C(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                    C(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                        +volume*sh(t1)/Mat_state.k(i);
                end
            end
        end
    end
end