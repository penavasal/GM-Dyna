classdef DYN_MATRIX
    
    properties 
        
        mass;
        damp;
        l_mass;
        l_damp;
        l_mass_w;
        l_mass_wn;
        
    end
    
    methods
        function obj=DYN_MATRIX
            
            global GEOMETRY SOLVER 


            sp=GEOMETRY.sp;
            df=GEOMETRY.sp;

            obj.l_mass=zeros(GEOMETRY.nodes*sp);
            obj.l_damp=zeros(GEOMETRY.nodes*sp);
            
            obj.damp=zeros(GEOMETRY.nodes*df);
            obj.mass=zeros(GEOMETRY.nodes*df);
            
            if SOLVER.UW
                obj.l_mass_w=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
                if SOLVER.IMPLICIT==0
                    obj.l_mass_wn=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
                end
            end
            
        end
    end
    
    methods (Static)  
        
        function [obj]=matrices(Mat_state,MAT_POINT,d,obj)

            global MATERIAL TIME SOLVER GEOMETRY VARIABLE

            alpha=TIME.alpha;
            rho_w=VARIABLE.rho_w;
            sp=GEOMETRY.sp;
            df=GEOMETRY.df;

            mass_mtx=zeros(GEOMETRY.nodes*df);
            damp_mtx=zeros(GEOMETRY.nodes*df);

            if alpha || SOLVER.UW==1  || SOLVER.UW==2
                for i=1:GEOMETRY.mat_points
                    
                    %Volume & density
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
                    
                    % Shape function
                    nd = MAT_POINT(i).near;
                    m  = length(nd);
                    sh = MAT_POINT(i).N;
                    if SOLVER.UW==2
                        b = MAT_POINT(i).B;
                        if SOLVER.AXI
                            mm=[1 1 0 1];
                        else
                            mm=[1 1 0];
                        end
                        Qt=(b'*mm'*sh')';

                        K_w=MATERIAL.MAT(28,MATERIAL.e(i));
                        K_s=MATERIAL.MAT(27,MATERIAL.e(i));

                        Q=1/(n/K_w+(1-n)/K_s);
                        
                        [A]=DYN_MATRIX.A_mat(m,sh,b);
                    end

                    for t1=1:m
                        for t2=1:m
                            for k=1:sp
                                if SOLVER.UW==1
                                    damp_mtx(nd(t1)*df+1-k,nd(t2)*df+1-k)=...
                                        damp_mtx(nd(t1)*df+1-k,nd(t2)*df+1-k)-...
                                        t*sh(t1)*sh(t2)/Mat_state.k(i);
                                elseif SOLVER.UW==2
                                	damp_mtx(nd(t1)*df,nd(t2)*df-k)=...
                                        damp_mtx(nd(t1)*df,nd(t2)*df-k)+...
                                        t*Qt(t1,t2*sp+1-k);
                                    damp_mtx(nd(t1)*df,nd(t2)*df)=...
                                        damp_mtx(nd(t1)*df,nd(t2)*df)+...
                                        t*sh(t1)*sh(t2)/Q;
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
                                elseif SOLVER.UW==2 && alpha
                                    mass_mtx(nd(t1)*df-k,nd(t2)*df-k)=...
                                        mass_mtx(nd(t1)*df-k,nd(t2)*df-k)-...
                                        t*sh(t1)*sh(t2)*dens;
                                	mass_mtx(nd(t1)*df,nd(t2)*df-k)=...
                                        mass_mtx(nd(t1)*df,nd(t2)*df-k)-...
                                        t*A(t1,t2*sp+1-k)*rho_w*Mat_state.k(i);
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
                    [M1]=DYN_MATRIX.mass_lin_uw(MAT_POINT,d);
                    mass_mtx=mass_mtx+M1;
                end
                [C1]=DYN_MATRIX.damp_lin_uw(MAT_POINT,Mat_state,d);
                damp_mtx=damp_mtx+C1;
            end
            
            obj.damp=damp_mtx;
            obj.mass=mass_mtx;

        end

        function [mass_mtx]=mass_lin_uw(MAT_POINT,d)

            global GEOMETRY VARIABLE MATERIAL SOLVER

            sp=GEOMETRY.sp;
            df=GEOMETRY.df;
            rho_w=VARIABLE.rho_w;

            mass_mtx=zeros(GEOMETRY.nodes*GEOMETRY.df);

            for i=1:GEOMETRY.mat_points

                volume=GEOMETRY.Area(i)*MAT_POINT(i).J;
                
                if SOLVER.AXI
                    t=2*pi*MAT_POINT(i).xg(1)*volume;
                else
                    t=volume;
                end
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

                MAT1=t*rho_w*N'*(U+W)*T;
                MAT2=t*rho_w*N'*U*T;
                MAT3=t*rho_w*(2*n-1)/n^2*N'*W*T;
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

            global GEOMETRY MATERIAL SOLVER

            sp=GEOMETRY.sp;
            df=GEOMETRY.df;

            damp_mtx=zeros(GEOMETRY.nodes*df);

            for i=1:GEOMETRY.mat_points

                volume=GEOMETRY.Area(i)*MAT_POINT(i).J;
                n=1-(1-MATERIAL.MAT(16,MATERIAL.e(i)))/MAT_POINT(i).J;
                
                if SOLVER.AXI
                    t=2*pi*MAT_POINT(i).xg(1)*volume;
                else
                    t=volume;
                end

                nd = MAT_POINT(i).near;
                m  = length(nd);
                sh = MAT_POINT(i).N;
                b  = MAT_POINT(i).B;

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

                MAT3=t/Mat_state.k(i)*N'*W*(1-(1-n)/perm(i)*dk)*T;        
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

        function [obj]=lumped_damp(MAT_POINT,Mat_state,obj)

            global GEOMETRY SOLVER

            sp=GEOMETRY.sp;

            C=zeros(GEOMETRY.nodes*sp);
            if SOLVER.UW
                % Lumped Damp **********************
                for i=1:GEOMETRY.mat_points
                    volume=GEOMETRY.Area(i)*MAT_POINT(i).J;
                    
                    if SOLVER.AXI
                        t=2*pi*MAT_POINT(i).xg(1)*volume;
                    else
                        t=volume;
                    end
                    nd = MAT_POINT(i).near;
                    m  = length(nd);
                    sh = MAT_POINT(i).N;
                    for t1=1:m
                        for k=1:sp
                            C(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                            C(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                                +t*sh(t1)/Mat_state.k(i);
                        end
                    end
                end
            end
            obj.l_damp=C;
        end
        
        function [obj]=lumped_mass(MAT_POINT,obj)

            global MATERIAL GEOMETRY SOLVER VARIABLE

            Material=MATERIAL.e;
            MAT=MATERIAL.MAT;

            sp=GEOMETRY.sp;


            mass=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);

            if SOLVER.UW
                mass_w=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
                if SOLVER.IMPLICIT==0
                    mass_wn=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
                end
            end

            % Lumped Mass **********************
            for i=1:GEOMETRY.mat_points
                volume=GEOMETRY.Area(i)*MAT_POINT(i).J;
                if SOLVER.UW
                    n=1-(1-MAT(16,Material(i)))/MAT_POINT(i).J;
                    dens=n*VARIABLE.rho_w+(1-n)*MAT(3,Material(i));
                else
                    dens=MAT(3,Material(i))/MAT_POINT(i).J;
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
                    for k=1:sp
                        mass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                        mass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                            +dens*t*sh(t1);

                        if SOLVER.UW
                            mass_w(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                            mass_w(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                                +VARIABLE.rho_w*t*sh(t1);
                            if SOLVER.IMPLICIT==0
                                mass_wn(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                                mass_wn(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                                    +VARIABLE.rho_w/n*t*sh(t1);
                            end
                        end
                    end
                end
            end
            
            obj.l_mass=mass;
            if SOLVER.UW
                obj.l_mass_w=mass_w;
                if SOLVER.IMPLICIT==0
                    obj.l_mass_wn=mass_wn;
                end
            end
            
        end
        
        function [obj]=lumped_mass_bf(MAT_POINT,Mat_state,obj)

            global MATERIAL GEOMETRY SOLVER VARIABLE

            Material=MATERIAL.e;
            MAT=MATERIAL.MAT;

            sp=GEOMETRY.sp;
            df=GEOMETRY.df;

            mass=zeros(GEOMETRY.nodes*df,GEOMETRY.nodes*df);

            % Lumped Mass **********************
            for i=1:GEOMETRY.mat_points
                volume=GEOMETRY.Area(i)*MAT_POINT(i).J;
                if SOLVER.UW
                    n=1-(1-MAT(16,Material(i)))/MAT_POINT(i).J;
                    dens=n*VARIABLE.rho_w+(1-n)*MAT(3,Material(i));
                else
                    dens=MAT(3,Material(i))/MAT_POINT(i).J;
                end
                
                if SOLVER.AXI
                    t=2*pi*MAT_POINT(i).xg(1)*volume;
                else
                    t=volume;
                end

                nd = MAT_POINT(i).near;
                m  = length(nd);
                sh = MAT_POINT(i).N;
                
                if SOLVER.UW==2
                    b = MAT_POINT(i).B;                  
                    dN=zeros(2,m);
                    for j=1:m
                        dN(1,j)=b(1,(j-1)*sp+1);
                        dN(2,j)=b(2,(j-1)*sp+2);
                    end
                end
                
 
                for t1=1:m
                    for k=1:sp
                        if SOLVER.UW==1
                            mass(nd(t1)*df-sp+1-k,nd(t1)*df-sp+1-k)=...
                            mass(nd(t1)*df-sp+1-k,nd(t1)*df-sp+1-k)...
                                +dens*t*sh(t1);
                            mass(nd(t1)*df+1-k,nd(t1)*df+1-k)=...
                            mass(nd(t1)*df+1-k,nd(t1)*df+1-k)...
                                +VARIABLE.rho_w*t*sh(t1);
                        elseif SOLVER.UW==2
                            mass(nd(t1)*df-k,nd(t1)*df-k)=...
                            mass(nd(t1)*df-k,nd(t1)*df-k)...
                                +dens*t*sh(t1);
                            mass(nd(t1)*df,nd(t1)*df-k)=...
                            mass(nd(t1)*df,nd(t1)*df-k)...
                                +VARIABLE.rho_w*t*dN(1+sp-k,t1)*Mat_state.k(i);
                        else
                            mass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                                mass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                                +dens*t*sh(t1);
                        end
                    end
                end
            end
            obj.l_mass=mass;           
        end
        
        function [A]=A_mat(n,N,b)
                   
            sp=2;
            
            dN=zeros(2,n);
            Nv=zeros(2,2*n);
            for j=1:n
                dN(1,j)=b(1,(j-1)*sp+1);
                dN(2,j)=b(2,(j-1)*sp+2);

                Nv(1,2*j-1)=N(j,1);
                Nv(2,2*j)=N(j,1);
            end

            A=dN'*Nv;         
        end
        
    end
end