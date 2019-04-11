function [mass_mtx,damp_mtx]=dyn_matrices(Mat_state,MAT_POINT,d)

    global MATERIAL TI_param SOLVER GEOMETRY VARIABLE
    
    alpha=TI_param(4);
    rho_w=VARIABLE.rho_w;
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    
    mass_mtx=zeros(GEOMETRY.nodes*df,GEOMETRY.nodes*df);
    damp_mtx=zeros(GEOMETRY.nodes*df,GEOMETRY.nodes*df);
    
    if alpha || SOLVER.UW==1 
        for i=1:GEOMETRY.mat_points
            volume=GEOMETRY.Area(i)*MAT_POINT(i).J;
            if SOLVER.UW
                n=1-(1-MATERIAL.MAT(16,MATERIAL.e(i)))/MAT_POINT(i).J;
                dens=n*rho_w+(1-n)*MATERIAL.MAT(3,MATERIAL.e(i));
            else
                dens=MATERIAL.MAT(3,MATERIAL.e(i))/MAT_POINT(i).J;
            end

            if SOLVER.AXI
                t=2*pi*MAT_POINT(i).xg(1)*volume;
            else
                t=volume;
            end
            nd = MAT_POINT(i).near;
            m  = length(nd);
            sh = MAT_POINT(i).N;

            for t1=1:m
                for t2=1:m
                    for k=1:sp
                        if SOLVER.UW==1
                            damp_mtx(nd(t1)*df+1-k,nd(t2)*df+1-k)=...
                                damp_mtx(nd(t1)*df+1-k,nd(t2)*df+1-k)-...
                                t*sh(t1)*sh(t2)/Mat_state.k(i);
                        end
                        if SOLVER.UW==1 && alpha
                            mass_mtx(nd(t1)*df-1-k,nd(t2)*df-1-k)=...
                                mass_mtx(nd(t1)*df-1-k,nd(t2)*df-1-k)-...
                                t*sh(t1)*sh(t2)*dens;
                            mass_mtx(nd(t1)*df-1-k,nd(t2)*df+1-k)=...
                                mass_mtx(nd(t1)*df-1-k,nd(t2)*df+1-k)-...
                                t*sh(t1)*sh(t2)*rho_w;
                            mass_mtx(nd(t1)*df+1-k,nd(t2)*df-1-k)=...
                                mass_mtx(nd(t1)*df+1-k,nd(t2)*df-1-k)-...
                                t*sh(t1)*sh(t2)*rho_w;
                            mass_mtx(nd(t1)*df+1-k,nd(t2)*df+1-k)=...
                                mass_mtx(nd(t1)*df+1-k,nd(t2)*df+1-k)-...
                                t*sh(t1)*sh(t2)*rho_w/n;
                        elseif alpha
                            mass_mtx(nd(t1)*df+1-k,nd(t2)*df+1-k)=...
                                mass_mtx(nd(t1)*df+1-k,nd(t2)*df+1-k)+...
                                +t*dens*sh(t1)*sh(t2);
                        end
                    end
                end
            end
        end
    end
    
    if SOLVER.LIN && SOLVER.UW==1 && SOLVER.AXI==0
        if alpha
            [M1]=mass_lin_uw(MAT_POINT,d);
            mass_mtx=mass_mtx+M1;
        end
        [C1]=damp_lin_uw(MAT_POINT,Mat_state,d);
        damp_mtx=damp_mtx+C1;
    end
    
end

function [mass_mtx]=mass_lin_uw(MAT_POINT,d)

    global GEOMETRY VARIABLE MATERIAL
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    rho_w=VARIABLE.rho_w;
    
    mass_mtx=zeros(GEOMETRY.nodes*GEOMETRY.df);
    
    for i=1:GEOMETRY.mat_points
        
        volume=GEOMETRY.Area(i)*MAT_POINT(i).J;
        n=1-(1-MATERIAL.MAT(16,MATERIAL.e(i)))/MAT_POINT(i).J;

        nd = MAT_POINT(i).near;
        m  = length(nd);
        sh = MAT_POINT(i).N;
        b  = MAT_POINT(i).B;
        
        T=[1 1 0]*b;
        
        u=zeros(m*sp,1);     % displacements _ solid
        w=zeros(m*sp,1);     % displacements _ water

        for j=1:m
            for l=1:sp
                u((i-1)*sp+l,1)=d((nd(i)-1)*df+l,1);
                w((i-1)*sp+l,1)=d((nd(i)-1)*df+sp+l,1);
            end
        end
        
        N=zeros(GEOMETRY.sp,m*GEOMETRY.sp);
        for j=1:m
            N(1,(j-1)*GEOMETRY.sp+1)=sh(j);
            N(2,(j-1)*GEOMETRY.sp+2)=sh(j);
        end
        
        U=N*u;
        W=N*w;
        
        MAT1=volume*rho_w*N'*(U+W)*T;
        MAT2=volume*rho_w*N'*U*T;
        MAT3=volume*rho_w*(2*n-1)/n^2*N'*W*T;
        MAT3=MAT2+MAT3;
        
        mat=zeros(m*df,m*df);
        for j=1:m
            for r=1:m
                for k=1:sp
                    for l=1:sp
                        mat((j-1)*df+k,(r-1)*df+l)=...
                            MAT1((j-1)*sp+k,(r-1)*sp+l);
                        mat((j-1)*df+k+2,(r-1)*df+l)=...
                            MAT3((j-1)*sp+k,(r-1)*sp+l);
                    end
                end
            end
        end

        % ----------------------------
        % Assembling of M
        % ----------------------------     
        for j=1:m
            for l=1:m
                for k=1:df
                    for r=1:df
                        mass_mtx(nd(j)*df+1-k,nd(l)*df+1-r)...
                        =mass_mtx(nd(j)*df+1-k,nd(l)*df+1-r)...
                        -mat(j*df+1-k,l*df+1-r);
                    end
                end
            end
        end
        
        clear MAT mat U W u w N
        
    end
end

function [damp_mtx]=damp_lin_uw(MAT_POINT,Mat_state,d)
     
    global GEOMETRY MATERIAL
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;

    damp_mtx=zeros(GEOMETRY.nodes*df);
    
    for i=1:GEOMETRY.mat_points
        
        volume=GEOMETRY.Area(i)*MAT_POINT(i).J;
        n=1-(1-MATERIAL.MAT(16,MATERIAL.e(i)))/MAT_POINT(i).J;

        nd = Shape_function.near{i};
        m  = length(nd);
        sh = Shape_function.p{i};
        b  = Shape_function.B{i};
        
        dk=0;
        
        T=[1 1 0]*b;
        w=zeros(m*sp,1);     % displacements _ water

        for j=1:m
            for l=1:sp
                w((i-1)*sp+l,1)=d((nd(i)-1)*df+sp+l,1);
            end
        end
        
        N=zeros(sp,m*sp);
        for j=1:m
            N(1,(j-1)*sp+1)=sh(j);
            N(2,(j-1)*sp+2)=sh(j);
        end
        
        W=N*w;
        
        MAT3=volume/Mat_state.k(i)*N'*W*(1-(1-n)/perm(i)*dk)*T;        
        mat=zeros(m*df,m*df);
        for j=1:m
            for r=1:m
                for k=1:sp
                    for l=1:sp
                        mat((j-1)*df+k+2,(r-1)*df+l)=...
                            MAT3((j-1)*sp+k,(r-1)*sp+l);
                    end
                end
            end
        end

        % ----------------------------
        % Assembling of C
        % ----------------------------     
        for j=1:m
            for l=1:m
                for k=1:df
                    for r=1:df
                        damp_mtx(nd(j)*df+1-k,nd(l)*df+1-r)...
                        =damp_mtx(nd(j)*df+1-k,nd(l)*df+1-r)...
                        -mat(j*df+1-k,l*df+1-r);
                    end
                end
            end
        end
        
        clear MAT mat W w N
        
    end
    
end