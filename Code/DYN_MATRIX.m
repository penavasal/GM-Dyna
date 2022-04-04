classdef DYN_MATRIX
    
    properties 
        
        mass;
        damp;
        l_mass;
        l_damp;
        l_damp_w;
        l_mass_w;
        l_mass_wn;
        l_mass_p;
        
    end
    
    methods
        function obj=DYN_MATRIX
            
            global GEOMETRY SOLVER 


            sp=GEOMETRY.sp;
            df=GEOMETRY.sp;

            obj.l_mass=zeros(GEOMETRY.nodes*sp);
            obj.l_damp=zeros(GEOMETRY.nodes*sp);
            obj.l_damp_w=zeros(GEOMETRY.nodes*sp);
            
            obj.damp=zeros(GEOMETRY.nodes*df);
            obj.mass=zeros(GEOMETRY.nodes*df);
            
            if SOLVER.UW
                obj.l_mass_w=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
                if SOLVER.IMPLICIT==0
                    obj.l_mass_wn=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
                    if SOLVER.UW==3
                        obj.l_mass_p=zeros(GEOMETRY.nodes,GEOMETRY.nodes);
                    end
                end
            end
            
        end
    end
    
    methods (Static)  
        
        function [obj]=matrices(Mat_state,MAT_POINT,d,obj,BLCK)

            global MATERIAL TIME SOLVER GEOMETRY BOUNDARY

            alpha=TIME{BLCK}.alpha;
            sp=GEOMETRY.sp;
            df=GEOMETRY.df;
            
            mati=GEOMETRY.material;

            mass_mtx=zeros(GEOMETRY.nodes*df);
            damp_mtx=zeros(GEOMETRY.nodes*df);
            
            
            
            if SOLVER.absorbing 
                constrains  = BOUNDARY{BLCK}.constrains;
                [dofs,bds]=size(constrains);
                abc=zeros(dofs,1);
                for i=1:dofs
                    for j=1:bds
                        if constrains(i,j)==4
                            abc(i)=1;
                            break;
                        end
                    end
                end
            end

            if (alpha || SOLVER.UW==1 || SOLVER.UW==4  || SOLVER.UW==2 || SOLVER.UW==3) && ...
                    SOLVER.DYN(BLCK)==1
                for i=1:GEOMETRY.mat_points
                    
                    rho_w=MATERIAL(BLCK).MAT{42,mati(i)};
                    
                    %Volume & density
                    volume=GEOMETRY.Area(i)*MAT_POINT{1}(i).J;
                    if SOLVER.UW
                        if SOLVER.UW==4
                            sw=Mat_state.sw(i,1);
                        else
                            sw=1;
                        end
                        n=1-(1-MATERIAL(BLCK).MAT{16,mati(i)})/MAT_POINT{1}(i).J;
                        dens=n*rho_w+sw*(1-n)*MATERIAL(BLCK).MAT{3,mati(i)};
                    else
                        dens=MATERIAL(BLCK).MAT{3,mati(i)}/MAT_POINT{1}(i).J;
                    end

                    if SOLVER.AXI
                        t=2*pi*MAT_POINT{1}(i).xg(1)*volume;
                    else
                        t=volume;
                    end
                     
                    
                    % Shape function
                    nd = MAT_POINT{1}(i).near;
                    m  = length(nd);
                    sh = MAT_POINT{1}(i).N;                
                    
                    if SOLVER.UW>0
                        ndw = MAT_POINT{2}(i).near;
                        mw  = length(ndw);
                        shw = MAT_POINT{2}(i).N;
                    end
                    
                    %Absorbing
                    if SOLVER.absorbing
                        for j=1:m
                            for k=1:df
                                if abc(nd(j)*df+1-k)==1
                                    
                                end
                            end
                        end
                    end
                    
                    if SOLVER.UW==3
                        ndp = MAT_POINT{3}(i).near;
                        mp  = length(ndp);
                        shp = MAT_POINT{3}(i).N;
                        b = MAT_POINT{1}(i).B;
                        if SOLVER.AXI
                            mm=[1 1 0 1];
                        else
                            mm=[1 1 0];
                        end
                        Qt=(b'*mm'*shp')';
                        
                        K_w=MATERIAL(BLCK).MAT{28,mati(i)};
                        K_s=MATERIAL(BLCK).MAT{27,mati(i)};

                        Q=1/(n/K_w+(1-n)/K_s);
                        
                        % stab
                        if SOLVER.Pstab
                            bp = MAT_POINT{1}(i).B;
                            dN=zeros(2,mp);
                            for j=1:mp
                                dN(1,j)=bp(1,(j-1)*sp+1);
                                dN(2,j)=bp(2,(j-1)*sp+2);
                            end
                            h=GEOMETRY.h_ini(i);
                            tau=SOLVER.Pstab*h*h/Q;
                            Hs=tau*(dN'*dN);
                        else
                            Hs=zeros(mp);
                        end  
                        
                    elseif SOLVER.UW==2
                        bw= MAT_POINT{2}(i).B;
                        b = MAT_POINT{1}(i).B;
                        if SOLVER.AXI
                            mm=[1 1 0 1];
                        else
                            mm=[1 1 0];
                        end
                        Qt=(b'*mm'*shw')';

                        K_w=MATERIAL(BLCK).MAT{28,mati(i)};
                        K_s=MATERIAL(BLCK).MAT{27,mati(i)};

                        Q=1/(n/K_w+(1-n)/K_s);
                        
                        [A]=DYN_MATRIX.A_mat(m,mw,sh,bw);
                        
                        % stab
                        if SOLVER.Pstab
                            dN=zeros(2,mw);
                            for j=1:mw
                                dN(1,j)=bw(1,(j-1)*sp+1);
                                dN(2,j)=bw(2,(j-1)*sp+2);
                            end
                            h=GEOMETRY.h_ini(i);
                            tau=SOLVER.Pstab*h*h/Q;
                            Hs=tau*(dN'*dN);
                        else
                            Hs=zeros(mw);
                        end   
                        
                    end

                    for t1=1:m
                        for t2=1:m
                            for k=1:sp
                                if (SOLVER.UW==1 || SOLVER.UW==4) && alpha
                                    mass_mtx(nd(t1)*df-1-k,nd(t2)*df-1-k)=...
                                        mass_mtx(nd(t1)*df-1-k,nd(t2)*df-1-k)-...
                                        t*sh(t1)*sh(t2)*dens;
                                elseif SOLVER.UW==2 && alpha
                                    mass_mtx(nd(t1)*df-k,nd(t2)*df-k)=...
                                        mass_mtx(nd(t1)*df-k,nd(t2)*df-k)-...
                                        t*sh(t1)*sh(t2)*dens;
                                elseif SOLVER.UW==3 && alpha
                                    mass_mtx(nd(t1)*df-sp-k,nd(t2)*df-sp-k)=...
                                        mass_mtx(nd(t1)*df-sp-k,nd(t2)*df-sp-k)-...
                                        t*sh(t1)*sh(t2)*dens;
                                elseif alpha
                                    mass_mtx(nd(t1)*df+1-k,nd(t2)*df+1-k)=...
                                        mass_mtx(nd(t1)*df+1-k,nd(t2)*df+1-k)+...
                                        +t*dens*sh(t1)*sh(t2);
                                end
                            end
                        end
                    end
                    if SOLVER.UW>0
                        for t1=1:mw
                            for t2=1:m
                                for k=1:sp
                                    if SOLVER.UW==2
                                        damp_mtx(ndw(t1)*df,nd(t2)*df-k)=...
                                            damp_mtx(ndw(t1)*df,nd(t2)*df-k)-...
                                            t*Qt(t1,t2*sp+1-k);
                                    end
                                    if (SOLVER.UW==1 || SOLVER.UW==4) && alpha
                                        mass_mtx(ndw(t1)*df-1-k,nd(t2)*df+1-k)=...
                                            mass_mtx(ndw(t1)*df-1-k,nd(t2)*df+1-k)-...
                                            t*shw(t1)*sh(t2)*rho_w;
                                    elseif SOLVER.UW==2 && alpha
                                        mass_mtx(ndw(t1)*df,nd(t2)*df-k)=...
                                            mass_mtx(ndw(t1)*df,nd(t2)*df-k)+...
                                            t*A(t1,t2*sp+1-k)*rho_w*Mat_state.k(i);
                                    elseif SOLVER.UW==3 && alpha
                                        mass_mtx(ndw(t1)*df-sp-k,nd(t2)*df-k)=...
                                            mass_mtx(ndw(t1)*df-sp-k,nd(t2)*df-k)-...
                                            t*shw(t1)*sh(t2)*rho_w;
                                    end
                                end
                            end
                        end
                        for t1=1:m
                            for t2=1:mw
                                for k=1:sp
                                    if (SOLVER.UW==1 || SOLVER.UW==4)  && alpha
                                        mass_mtx(nd(t1)*df+1-k,ndw(t2)*df-1-k)=...
                                            mass_mtx(nd(t1)*df+1-k,ndw(t2)*df-1-k)-...
                                            t*sh(t1)*shw(t2)*rho_w;
                                    elseif SOLVER.UW==3  && alpha
                                        mass_mtx(nd(t1)*df-k,ndw(t2)*df-sp-k)=...
                                            mass_mtx(nd(t1)*df-k,ndw(t2)*df-sp-k)-...
                                            t*sh(t1)*shw(t2)*rho_w;
                                    end
                                end
                            end
                        end
                        for t1=1:mw
                            for t2=1:mw
                                for k=1:sp
                                    if SOLVER.UW==1 || SOLVER.UW==4
                                        damp_mtx(ndw(t1)*df+1-k,ndw(t2)*df+1-k)=...
                                            damp_mtx(ndw(t1)*df+1-k,ndw(t2)*df+1-k)-...
                                            t*shw(t1)*shw(t2)/Mat_state.k(i);
                                    elseif SOLVER.UW==3
                                        damp_mtx(ndw(t1)*df-k,ndw(t2)*df-k)=...
                                            damp_mtx(ndw(t1)*df-k,ndw(t2)*df-k)-...
                                            t*shw(t1)*shw(t2)/Mat_state.k(i);
                                    elseif SOLVER.UW==2
                                        damp_mtx(ndw(t1)*df,ndw(t2)*df)=...
                                            damp_mtx(ndw(t1)*df,ndw(t2)*df)-...
                                            t*shw(t1)*shw(t2)/Q+t*Hs(t1,t2);
                                    end
                                    if (SOLVER.UW==1 || SOLVER.UW==4) && alpha
                                        mass_mtx(ndw(t1)*df+1-k,ndw(t2)*df+1-k)=...
                                            mass_mtx(ndw(t1)*df+1-k,ndw(t2)*df+1-k)-...
                                            t*shw(t1)*shw(t2)*rho_w/n/sw;
                                    elseif SOLVER.UW==3 && alpha
                                        mass_mtx(ndw(t1)*df-k,ndw(t2)*df-k)=...
                                            mass_mtx(ndw(t1)*df-k,ndw(t2)*df-k)-...
                                            t*shw(t1)*shw(t2)*rho_w/n/sw;
                                    end
                                end
                            end
                        end
                        if SOLVER.UW==3
                            for t1=1:mp
                                for t2=1:m
                                    for k=1:sp
                                        damp_mtx(ndp(t1)*df,nd(t2)*df-sp-k)=...
                                            damp_mtx(ndp(t1)*df,nd(t2)*df-sp-k)-...
                                            t*Qt(t1,t2*sp+1-k);
                                    end
                                end
                                for t2=1:mw
                                    for k=1:sp
                                        damp_mtx(ndp(t1)*df,ndw(t2)*df-k)=...
                                            damp_mtx(ndp(t1)*df,ndw(t2)*df-k)-...
                                            t*Qt(t1,t2*sp+1-k);
                                    end
                                end
                                for t2=1:mp
                                    damp_mtx(ndw(t1)*df,ndw(t2)*df)=...
                                        damp_mtx(ndw(t1)*df,ndw(t2)*df)-...
                                        t*shw(t1)*shw(t2)/Q-t*Hs(t1,t2);
                                end
                            end
                        end
                    end
                end
            end

            if SOLVER.LIN && (SOLVER.UW==1 || SOLVER.UW==4) && SOLVER.AXI==0 && SOLVER.DYN(BLCK)==1
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

        function [mass_mtx]=mass_lin_uw(MAT_POINT,d,BLCK)

            global GEOMETRY VARIABLE MATERIAL SOLVER

            sp=GEOMETRY.sp;
            df=GEOMETRY.df;
            rho_w=VARIABLE.rho_w;
            
            mati=GEOMETRY.material;

            mass_mtx=zeros(GEOMETRY.nodes*GEOMETRY.df);

            for i=1:GEOMETRY.mat_points

                volume=GEOMETRY.Area(i)*MAT_POINT{1}(i).J;
                
                if SOLVER.AXI
                    t=2*pi*MAT_POINT{1}(i).xg(1)*volume;
                else
                    t=volume;
                end
                n=1-(1-MATERIAL(BLCK).MAT(16,mati(i)))/MAT_POINT{1}(i).J;

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

        function [damp_mtx]=damp_lin_uw(MAT_POINT,Mat_state,d,BLCK)

            global GEOMETRY MATERIAL SOLVER

            sp=GEOMETRY.sp;
            df=GEOMETRY.df;
            
            mati=GEOMETRY.material;

            damp_mtx=zeros(GEOMETRY.nodes*df);

            for i=1:GEOMETRY.mat_points

                volume=GEOMETRY.Area(i)*MAT_POINT{1}(i).J;
                n=1-(1-MATERIAL(BLCK).MAT(16,mati(i)))/MAT_POINT{1}(i).J;
                
                if SOLVER.AXI
                    t=2*pi*MAT_POINT{1}(i).xg(1)*volume;
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

        function [obj]=expl_mat(MAT_POINT,Mat_state,obj,BLCK,STEP)
            
            global MATERIAL GEOMETRY SOLVER BOUNDARY

            Material=GEOMETRY.material;
            MAT=MATERIAL(BLCK).MAT;

            sp=GEOMETRY.sp;
            t   = STEP.t;

            % Allocate
            mass=zeros(GEOMETRY.nodes*sp);
            if SOLVER.UW==1||SOLVER.UW==4
                Cw=zeros(GEOMETRY.nodes*sp);
                mass_w=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
                if SOLVER.IMPLICIT==0 
                    mass_wn=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
                end
            elseif SOLVER.UW==2
                Cw=zeros(GEOMETRY.nodes,GEOMETRY.nodes*sp);
                mass_w=zeros(GEOMETRY.nodes,GEOMETRY.nodes*sp);
            elseif SOLVER.UW==3
                Cw=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
                mass_w=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
                mass_p=zeros(GEOMETRY.nodes,GEOMETRY.nodes);
                mass_wn=zeros(GEOMETRY.nodes*sp,GEOMETRY.nodes*sp);
            end

            % Absorbing BC
            if SOLVER.absorbing
                C=zeros(GEOMETRY.nodes*sp);
                ABC    = BOUNDARY{BLCK}.abc;
                bds    = BOUNDARY{BLCK}.size;
                type   = BOUNDARY{BLCK}.Type;
                b_mult = BOUNDARY{BLCK}.b_mult; 
                
                LIST=[];
                for j=1:bds
                    if t<eval(b_mult(2,j)) || t>eval(b_mult(3,j))
                        continue;
                    elseif type(j)~=7
                        continue;
                    else
                        LIST= [LIST;ABC{j}];
                    end
                end
                [abcs,sz]=size(LIST);
            end
           
            
            
            for i=1:GEOMETRY.mat_points 
                % Volume
                volume=GEOMETRY.Area(i)*MAT_POINT{1}(i).J;
                if SOLVER.AXI
                    t=2*pi*MAT_POINT{1}(i).xg(1)*volume;
                else
                    t=volume;
                end
                
                % Density
                if SOLVER.UW
                    if SOLVER.UW==4 || SOLVER.UW==5
                        sw=Mat_state.sw(i,1);
                    else
                        sw=1;
                    end
                    rho_w=MAT{42,Material(i)};
                    n=1-(1-MAT{16,Material(i)})/MAT_POINT{1}(i).J;
                    dens=sw*n*rho_w+(1-n)*MAT{3,Material(i)};

                    if SOLVER.UW==3 || SOLVER.UW==2
                        K_w=MATERIAL(BLCK).MAT{28,Material(i)};
                        K_s=MATERIAL(BLCK).MAT{27,Material(i)};
                        
                        if sw==1
                            Q=1/(n/K_w+(1-n)/K_s);
                        else
                            Q=W_retention.Q(i,Mat_state,K_s,K_w,n);
                        end
                    end

                else
                    dens=MAT{3,Material(i)}/MAT_POINT{1}(i).J;
                end
                
                % Shape function
                nd = MAT_POINT{1}(i).near;
                m  = length(nd);
                sh = MAT_POINT{1}(i).N;
                
                % Absorbing
                if SOLVER.absorbing
                    elem=GEOMETRY.elem(MAT_POINT{1}(i).element,:);
                    j=1;
                    while j<=abcs
                        A=ismember(LIST(j,:),elem);
                        if A(1) && A(2)
                            V=zeros(sp,1);
                            for k=1:sp
                                V(k)=GEOMETRY.x_0(LIST(j,1),k)-...
                                    GEOMETRY.x_0(LIST(j,2),k);
                            end
                            d=sqrt(sum(V.^2));
                            nv=norm(V);
                            if nv==0
                                continue
                            else
                                V=V/norm(V);
                            end
                            if V(2)==0
                                theta=0;
                            elseif V(1)==0
                                theta=pi/2;
                            else
                                theta=atan(abs(V(2))/abs(V(1)));
                            end
                            vc=MAT{6,Material(i)};    
                            vs=sqrt(MAT{4,Material(i)}/dens); 
                            quad=LIST(j,sz);
                            if quad
                                d=d/4;
                            else
                                d=d/4;
                            end
                            
                            if SOLVER.UW==1 ||SOLVER.UW==4
                                K_w=MATERIAL(BLCK).MAT{28,Material(i)};
                                K_s=MATERIAL(BLCK).MAT{27,Material(i)};

                                Q=1/(n/K_w+(1-n)/K_s);
                                vcw=sqrt(Q/rho_w);
                            end
                            break;
                        end
                        j=j+1;
                    end
                    
                end
                
                % Solid Lumped Mass **********************
                for k=1:sp
                    for t1=1:m
                        mass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                        mass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                            +dens*t*sh(t1);
                    end
                end
                
                % Absorbing
                if SOLVER.absorbing && j<=abcs
                    for t1=1:m
                        for t2=1:sz
                            nda=LIST(j,t2);
                            if nd(t1)==nda
                                if quad
                                    if t2<sz
                                        C(nd(t1)*sp-1,nd(t1)*sp-1)=...
                                        C(nd(t1)*sp-1,nd(t1)*sp-1)...
                                            +dens*d*(vs*cos(theta)+vc*sin(theta));
                                        C(nd(t1)*sp,nd(t1)*sp)=C(nd(t1)*sp,nd(t1)*sp)...
                                            +dens*d*(vs*cos(pi/2-theta)+vc*sin(pi/2-theta));
                                    else
                                        C(nd(t1)*sp-1,nd(t1)*sp-1)=...
                                        C(nd(t1)*sp-1,nd(t1)*sp-1)...
                                            +2*dens*d*(vs*cos(theta)+vc*sin(theta));
                                        C(nd(t1)*sp,nd(t1)*sp)=C(nd(t1)*sp,nd(t1)*sp)...
                                            +2*dens*d*(vs*cos(pi/2-theta)+vc*sin(pi/2-theta));
                                    end
                                else
                                    C(nd(t1)*sp-1,nd(t1)*sp-1)=...
                                    C(nd(t1)*sp-1,nd(t1)*sp-1)...
                                        +dens*d*(vs*cos(theta)+vc*sin(theta));
                                    C(nd(t1)*sp,nd(t1)*sp)=C(nd(t1)*sp,nd(t1)*sp)...
                                        +dens*d*(vs*cos(pi/2-theta)+vc*sin(pi/2-theta));
                                end
                            end
                        end
                    end
                end
                
                if SOLVER.UW>0
                    ndw = MAT_POINT{2}(i).near;
                    mw  = length(ndw);
                    shw = MAT_POINT{2}(i).N;
                    if SOLVER.UW==3
                        ndp = MAT_POINT{3}(i).near;
                        mp  = length(ndp);
                        shp = MAT_POINT{3}(i).N;
                    end
                end
                
                if SOLVER.UW==1 || SOLVER.UW==4 || SOLVER.UW==3
                    val=1/Mat_state.k(i);
                    for t1=1:mw
                        for k=1:sp
                            % Lumped Damp **********************
                            Cw(ndw(t1)*sp+1-k,ndw(t1)*sp+1-k)=...
                            Cw(ndw(t1)*sp+1-k,ndw(t1)*sp+1-k)...
                                +t*shw(t1)*val;
                            % Lumped water mass **********************
                            mass_w(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                            mass_w(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                                +rho_w*t*shw(t1);
                            % Lumped water-porosity mass **********************
                            if SOLVER.IMPLICIT==0
                                mass_wn(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                                mass_wn(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                                    +rho_w/sw/n*t*shw(t1);
                            end
                        end
                        
                        % Absorbing water
                        if SOLVER.absorbing && j<=abcs
                            for t2=1:sz
                                nda=LIST(j,t2);
                                if ndw(t1)==nda
                                    if quad
                                        if t2<sz
                                            Cw(ndw(t1)*sp-1,ndw(t1)*sp-1)=...
                                            Cw(ndw(t1)*sp-1,ndw(t1)*sp-1)...
                                                +rho_w*d*vc*sin(theta);
                                            Cw(ndw(t1)*sp,ndw(t1)*sp)=Cw(ndw(t1)*sp,ndw(t1)*sp)...
                                                +rho_w*d*vc*sin(pi/2-theta);
                                        else
                                            Cw(ndw(t1)*sp-1,ndw(t1)*sp-1)=...
                                            Cw(ndw(t1)*sp-1,ndw(t1)*sp-1)...
                                                +2*rho_w*d*vc*sin(theta);
                                            Cw(ndw(t1)*sp,ndw(t1)*sp)=Cw(ndw(t1)*sp,ndw(t1)*sp)...
                                                +2*rho_w*d*vc*sin(pi/2-theta);
                                        end
                                    else
                                        Cw(ndw(t1)*sp-1,ndw(t1)*sp-1)=...
                                        Cw(ndw(t1)*sp-1,ndw(t1)*sp-1)...
                                            +rho_w*d*vc*sin(theta);
                                        Cw(ndw(t1)*sp,ndw(t1)*sp)=Cw(ndw(t1)*sp,ndw(t1)*sp)...
                                            +rho_w*d*vcw*sin(pi/2-theta);
                                    end
                                end
                            end
                        end
                    end

                    if SOLVER.UW==3
                        for t1=1:mp
                            % Lumped Damp **********************
                            mass_p(ndp(t1),ndp(t1))=...
                            mass_p(ndp(t1),ndp(t1))+t*shp(t1)/Q;
                        end

                    end

                elseif SOLVER.UW==2
                    bw= MAT_POINT{2}(i).B;
                    b = MAT_POINT{1}(i).B;
                    if SOLVER.AXI
                        mm=[1 1 0 1];
                    else
                        mm=[1 1 0];
                    end
                    Qt=(b'*mm'*shw')';

                    [A]=DYN_MATRIX.A_mat(m,mw,sh,bw);
                    
                    % stab
                    if SOLVER.Pstab
                        h=GEOMETRY.h_ini(i);
                        Vc=MATERIAL(BLCK).MAT{6,Material(i)};
                        tau=SOLVER.Pstab*h/Vc;
                    else
                        tau=0;
                    end 
                    
                    for t1=1:mw
                        for t2=1:m
                            for k=1:sp
                                Cw(ndw(t1),nd(t2)*sp+1-k)=...
                                    Cw(ndw(t1),nd(t2)*sp+1-k)-...
                                    t*Q*Qt(t1,t2*sp+1-k);

                                mass_w(ndw(t1),nd(t2)*sp+1-k)=...
                                    mass_w(ndw(t1),nd(t2)*sp+1-k)+...
                                    t*Q*A(t1,t2*sp+1-k)*(rho_w*Mat_state.k(i)-tau);
                            end
                        end
                    end
                end
            end
            
            obj.l_mass=mass;
            if SOLVER.absorbing
                obj.l_damp=C;
            end
            if SOLVER.UW
                obj.l_damp_w=Cw;
                obj.l_mass_w=mass_w;
                if SOLVER.IMPLICIT(BLCK)==0 && SOLVER.UW==1||SOLVER.UW==4||SOLVER.UW==3
                    obj.l_mass_wn=mass_wn;
                end
                if SOLVER.IMPLICIT(BLCK)==0 && SOLVER.UW==3
                    obj.l_mass_p=mass_p;
                end
            end
        end
        
        function [obj]=lumped_mass_bf(MAT_POINT,Mat_state,obj,BLCK)

            global MATERIAL GEOMETRY SOLVER

            Material=GEOMETRY.material;
            MAT=MATERIAL(BLCK).MAT;

            sp=GEOMETRY.sp;
            df=GEOMETRY.df;

            mass=zeros(GEOMETRY.nodes*df,GEOMETRY.nodes*df);

            % Lumped Mass **********************
            for i=1:GEOMETRY.mat_points
                volume=GEOMETRY.Area(i)*MAT_POINT{1}(i).J;
                if SOLVER.UW
                    if SOLVER.UW==4
                        sw=Mat_state.sw(i,1);
                    else
                        sw=1;
                    end
                    rho_w=MAT{42,Material(i)};
                    n=1-(1-MAT{16,Material(i)})/MAT_POINT{1}(i).J;
                    dens=n*rho_w+sw*(1-n)*MAT{3,Material(i)};
                else
                    dens=MAT{3,Material(i)}/MAT_POINT{1}(i).J;
                end
                
                if SOLVER.AXI
                    t=2*pi*MAT_POINT{1}(i).xg(1)*volume;
                else
                    t=volume;
                end

                nd = MAT_POINT{1}(i).near;
                m  = length(nd);
                sh = MAT_POINT{1}(i).N;
                
                if SOLVER.UW==1 || SOLVER.UW==4
                    ndw = MAT_POINT{2}(i).near;
                    mw  = length(ndw);
                    shw = MAT_POINT{2}(i).N;
                end
                
                if SOLVER.UW==2
                    b = MAT_POINT{2}(i).B;
                    ndw = MAT_POINT{2}(i).near;
                    mw  = length(ndw);
                    dN=zeros(2,mw);
                    shw = MAT_POINT{2}(i).N;
                    for j=1:mw
                        dN(1,j)=b(1,(j-1)*sp+1);
                        dN(2,j)=b(2,(j-1)*sp+2);
                    end
                end
                
                for k=1:sp
                    for t1=1:m
                        if SOLVER.UW==1 || SOLVER.UW==4
                            mass(nd(t1)*df-sp+1-k,nd(t1)*df-sp+1-k)=...
                            mass(nd(t1)*df-sp+1-k,nd(t1)*df-sp+1-k)...
                                +dens*t*sh(t1);
                        elseif SOLVER.UW==3
                            mass(nd(t1)*df-sp-k,nd(t1)*df-sp-k)=...
                            mass(nd(t1)*df-sp-k,nd(t1)*df-sp-k)...
                                +dens*t*sh(t1);
                        elseif SOLVER.UW==2
                            mass(nd(t1)*df-k,nd(t1)*df-k)=...
                            mass(nd(t1)*df-k,nd(t1)*df-k)...
                                +dens*t*sh(t1);
                        else
                            mass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                                mass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                                +dens*t*sh(t1);
                        end
                    end
                    if SOLVER.UW
                        for t1=1:mw
                            if SOLVER.UW==1 || SOLVER.UW==4
                                mass(ndw(t1)*df+1-k,ndw(t1)*df+1-k)=...
                                mass(ndw(t1)*df+1-k,ndw(t1)*df+1-k)...
                                    +rho_w*t*shw(t1);
                            elseif SOLVER.UW==3
                                mass(ndw(t1)*df-k,ndw(t1)*df-k)=...
                                mass(ndw(t1)*df-k,ndw(t1)*df-k)...
                                    +rho_w*t*shw(t1);
                            elseif SOLVER.UW==2
                                mass(ndw(t1)*df,ndw(t1)*df-k)=...
                                mass(ndw(t1)*df,ndw(t1)*df-k)...
                                    +rho_w*t*dN(1+sp-k,t1)*Mat_state.k(i);
                            end
                        end
                    end
                end
            end
            obj.l_mass=mass;           
        end
        
        function [A]=A_mat(n,nw,N,b)
                   
            sp=2;
            
            dN=zeros(2,n);
            Nv=zeros(2,2*n);
            for j=1:nw
                dN(1,j)=b(1,(j-1)*sp+1);
                dN(2,j)=b(2,(j-1)*sp+2);
            end
            for j=1:n
                Nv(1,2*j-1)=N(j,1);
                Nv(2,2*j)=N(j,1);
            end

            A=dN'*Nv;         
        end
        
    end
end
