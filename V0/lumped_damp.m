function [C]=lumped_damp(MAT_POINT,Mat_state)

    global GEOMETRY SOLVER
    
    sp=GEOMETRY.sp;

    C=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
    if SOLVER.UW
        %% Lumped Damp **********************
        for i=1:GEOMETRY.mat_points
            volume=GEOMETRY.Area(i)*MAT_POINT(i).J;
            nd = MAT_POINT(i).near;
            m  = length(nd);
            sh = MAT_POINT(i).N;
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