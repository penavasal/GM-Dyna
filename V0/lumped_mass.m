function [lumped]=lumped_mass(Shape_function,Mat_state)

    global MATERIAL GEOMETRY SOLVER VARIABLE
    
    Material=MATERIAL.e;
    MAT=MATERIAL.MAT;
    
    sp=GEOMETRY.sp;
    
    
    lumped.mass=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
    
    if SOLVER.UW
        lumped.mass_w=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
        if SOLVER.IMPLICIT==0
            lumped.mass_wn=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
        end
    end
    
    %% Lumped Mass **********************
    for i=1:GEOMETRY.elements
        volume=GEOMETRY.Area(i)*Mat_state.J(i);
        if SOLVER.UW
            n=1-(1-MAT(16,Material(i)))/Mat_state.J(i);
            dens=n*VARIABLE.rho_w+(1-n)*MAT(3,Material(i));
        else
            dens=MAT(3,Material(i))/Mat_state.J(i);
        end
        
        nd = Shape_function.near{i};
        m  = length(nd);
        sh = Shape_function.p{i};
        for t1=1:m
            for k=1:sp
                lumped.mass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                lumped.mass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                    +dens*volume*sh(t1);
                
                if SOLVER.UW
                    lumped.mass_w(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                    lumped.mass_w(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                        +VARIABLE.rho_w*volume*sh(t1);
                    if SOLVER.IMPLICIT==0
                        lumped.mass_wn(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                        lumped.mass_wn(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                            +VARIABLE.rho_w/n*volume*sh(t1);
                    end
                end
            end
        end
    end
end