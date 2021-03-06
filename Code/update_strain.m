
function [Mat_state,MAT_POINT]=update_strain(d,Mat_state,MAT_POINT,STEP)
    
    global SOLVER
    
    if SOLVER.SMALL==0
        [Mat_state,MAT_POINT]=update_F(d,Mat_state,MAT_POINT,STEP);
    else
        [Mat_state,MAT_POINT]=update_E(d,Mat_state,MAT_POINT);
    end
    
end

function [Mat_state,MAT_POINT]=update_F(d,Mat_state,MAT_POINT,STEP)
    
    global GEOMETRY SOLVER MATERIAL VARIABLE
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    dimf=GEOMETRY.f_dim;
    
    [aa,bb]=size(GEOMETRY.patch_con);
    area_p_n=zeros(aa,1);
    area_p_w=zeros(aa,1);
    logJ=zeros(aa,1);
    logJ_w=zeros(aa,1);

        
    m_I=LIB.m2v(eye(3));
    f_v=zeros(dimf,1);
    f_v_w=zeros(dimf,1);

    jacobians=zeros(GEOMETRY.mat_points,1);
        
    % Update **********************
    for j=1:aa
        for k=1:bb
            e=GEOMETRY.patch_con(j,k);

            % SOLID PHASE
            nd = MAT_POINT{1}(e).near;
            nn = length(nd);         
            b  = MAT_POINT{1}(e).B;
            if SOLVER.AXI
                bs=zeros(5,nn*2);
                for i=2:2:nn*2
                    bs(1,i-1)=b(1,i-1);
                    bs(2,i-1)=b(2,i);
                    bs(3,i)=b(3,i);
                    bs(4,i)=b(3,i-1);
                    bs(5,i-1)=b(4,i-1);
                end
                
            else
                bs=zeros(4,nn*2);
                for i=2:2:nn*2
                    bs(1,i-1)=b(1,i-1);
                    bs(2,i-1)=b(2,i);
                    bs(3,i)=b(3,i);
                    bs(4,i)=b(3,i-1);
                end
            end
            
            us=zeros(nn*sp,1);     %Incremental displacements _ solid
            for i=1:nn
                for l=1:sp
                    us((i-1)*sp+l,1)=d((nd(i)-1)*df+l,1)-d((nd(i)-1)*df+l,2);
                end
            end
            
            dF_=bs*us;                % Incremental F

            if SOLVER.AXI==0
                %dF_(dimf)=dF_(1);
                dF_(dimf)=0;
            end
            dF_=dF_+m_I;

            %Vector F to Matrix F
            for i=1:5
                f_v(i,1)=Mat_state.F((e-1)*5 + i,2);
            end           
            [F]=LIB.v2m(f_v);            
            [dF]=LIB.v2m(dF_);
            F=dF*F;           
            jacobians(e)=det(F);
            
            
            if SOLVER.UW==1 || SOLVER.UW==4
                ndw = MAT_POINT{2}(e).near;
                nn2 = length(ndw);         
                b2  = MAT_POINT{2}(e).B;
                if SOLVER.AXI
                    bw=zeros(5,nn2*2);
                    for i=2:2:nn2*2
                        bw(1,i-1)=b2(1,i-1);
                        bw(4,i)=b2(3,i-1);
                        bw(5,i-1)=b2(4,i-1);
                    end

                else
                    bw=zeros(4,nn2*2);
                    for i=2:2:nn2*2
                        bw(1,i-1)=b2(1,i-1);
                        bw(4,i)=b2(3,i-1);
                    end
                end
                uw=zeros(nn2*sp,1);     %Incremental displacements _ water
                for i=1:nn2
                    for l=1:sp
                        uw((i-1)*sp+l,1)=d((ndw(i)-1)*df+sp+l,1)...
                            -d((ndw(i)-1)*df+sp+l,2);
                    end
                end

                dF_w=bw*uw;
                if SOLVER.AXI==0
                    dF_w(dimf)=0;
                end
                dF_w=dF_w+m_I;

                %Vector F to Matrix F
                for i=1:dimf
                    f_v_w(i,1)=Mat_state.Fw((e-1)*dimf + i,2);
                end           
                [F_w]=LIB.v2m(f_v_w);
                [dFw]=LIB.v2m(dF_w);
                F_w=dFw*F_w;
                
            elseif SOLVER.UW==2 || SOLVER.UW==3
                if SOLVER.UW==2
                    N=MAT_POINT{2}(e).N;
                    ndw = MAT_POINT{2}(e).near;
                    nn2 = length(ndw);         
                    b2  = MAT_POINT{2}(e).B;
                else
                    N=MAT_POINT{3}(e).N;
                    ndw = MAT_POINT{3}(e).near;
                    nn2 = length(ndw);         
                    b2  = MAT_POINT{3}(e).B;
                end
                pwv=zeros(nn2,1);     %Incremental displacements _ solid
                dN=zeros(2,nn2);
                for i=1:nn2
                    dN(1,i)=b2(1,(i-1)*sp+1);
                    dN(2,i)=b2(2,(i-1)*sp+2);
                    pwv(i,1)=d(ndw(i)*df,1);
                end
                %Pw=N'*pwv;
                Mat_state.pw(e,1)=Mat_state.pw(e,3)+N'*pwv;
                dPw=dN*pwv;
                Mat_state.dpw((e-1)*sp+1:e*sp,1)=dPw;
            end
            
            if SOLVER.F_BAR>0
                %%%%%%   F_BAR   %%%%%%
                F_store(k)={F};
                % New vol, dens, jaco
                area_p_n(j)=area_p_n(j)+GEOMETRY.Area(e);%*jacobians(e);
                logJ(j)=logJ(j) + GEOMETRY.Area(e)*log(jacobians(e));%*jacobians(e);
            else
                % New jaco
                MAT_POINT{1}(e).J=jacobians(e);
                
                %Storage of vector F
                [f]=LIB.m2v(F);
                for i=1:dimf
                    Mat_state.F(e*dimf+1-i,1)=f(dimf+1-i);
                end

                if isnan(F)
                    disp('Fail in F');
                end
                
                if SOLVER.REMAPPING
                    %Principal strecht
                    C=F'*F;
                    [EP_]=Principal(C);
                    MAT_POINT{1}(e).EP(:,1)=EP_(:);
                end
            end
            if SOLVER.UW==1 || SOLVER.UW==4
                if SOLVER.F_BAR_W>0
                    %%%%%%   F_BAR W  %%%%%%
                    F_store_w(k)={F_w};
                    % New vol, dens, jaco
                    area_p_w(j)=area_p_w(j)+GEOMETRY.Area(e);%*jacobians(e);
                    logJ_w(j)=logJ_w(j) + GEOMETRY.Area(e)*log(det(F_w));%*jacobians(e);
                else
                    %Storage of vector F
                    [f_w]=LIB.m2v(F_w);
                    for i=1:dimf
                        Mat_state.Fw(e*dimf+1-i,1)=f_w(dimf+1-i);
                    end
                end
%             elseif SOLVER.UW==2
%                 if SOLVER.F_BAR_PW>0
%                     %%%%%%   F_BAR PW  %%%%%%
%                     F_store_w(k)={F_w};
%                     % New vol, dens, jaco
%                     area_p_w(j)=area_p_w(j)+GEOMETRY.Area(e)*det(F_w);
%                     logJ_w(j)=logJ_w(j) + GEOMETRY.Area(e)*det(F_w)*log(det(F_w));
%                 else
%                     %Storage of vector F
%                     [f_w]=LIB.m2v(F_w);
%                     for i=1:dimf
%                         Mat_state.Fw(e*dimf+1-i,1)=f_w(dimf+1-i);
%                     end
%                 end
            end
        end        
        %%%%%%   F_BAR   %%%%%%
        if SOLVER.F_BAR>0
            for k=1:bb
                F=F_store{k};
                e=GEOMETRY.patch_con(j,k);
                                
                J_=exp(SOLVER.F_BAR/area_p_n(j)*logJ(j)+...
                    (1-SOLVER.F_BAR)*log(jacobians(e)));
                
                if SOLVER.AXI
                    F_=(J_/jacobians(e))^(1/3)*F;
                else
                    F_=(J_/jacobians(e))^(1/2)*F(1:2,1:2);
                    F_(3,3)=1;
                end
                
                % New vol, dens, jaco
                jaco_new=det(F_);
                MAT_POINT{1}(e).J=jaco_new;
                
                %Storage of vector F
                [f]=LIB.m2v(F_);
                for i=1:dimf
                    Mat_state.F(e*dimf+1-i,1)=f(dimf+1-i);
                end
                
                if SOLVER.REMAPPING
                    %Principal strecht
                    C=F_'*F_;
                    [EP_]=Principal(C);
                    MAT_POINT{1}(e).EP(:,1)=EP_(:);
                end
            end
        end
        if SOLVER.UW==1 || SOLVER.UW==4
            if SOLVER.F_BAR_W>0
                for k=1:bb
                    F_w=F_store_w{k};
                    e=GEOMETRY.patch_con(j,k);

                    J_=exp(SOLVER.F_BAR_W/area_p_w(j)*logJ_w(j)+...
                        (1-SOLVER.F_BAR_W)*log(det(F_w)));


                    if SOLVER.AXI
                        F_w_=(J_/det(F_w))^(1/3)*F_w;
                    else
                        F_w_=(J_/det(F_w))^(1/2)*F_w(1:2,1:2);
                        F_w_(3,3)=1;
                    end

                    %F_w_=(J_bar/det(F_w))^(1/3)*F_w;

                    %Storage of vector F
                    [f_w]=LIB.m2v(F_w_);
                    for i=1:dimf
                        Mat_state.Fw(e*dimf+1-i,1)=f_w(dimf+1-i);
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        if SOLVER.UW>0
            mmat=MATERIAL(STEP.BLCK).MAT;
            for k=1:bb
                e=GEOMETRY.patch_con(j,k);
                mati=GEOMETRY.material(e);
                Mat_state.k(e) = ...
                    mmat{15,mati}/...
                    mmat{42,mati}/VARIABLE.g;
            end
        end
    end
end



function [Mat_state,MAT_POINT]=update_E(d,Mat_state,MAT_POINT)
    
    global GEOMETRY SOLVER
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    dims=GEOMETRY.s_dim;
    
        
    % Update **********************
    for e=1:GEOMETRY.mat_points

        % SOLID PHASE
        nd = MAT_POINT{1}(e).near;
        nn = length(nd);         
        b  = MAT_POINT{1}(e).B;
        
        us=zeros(nn*sp,1);     %Incremental displacements _ solid
        for i=1:nn
            for l=1:sp
                us((i-1)*sp+l,1)=d((nd(i)-1)*df+l,1);
            end
        end

        E_=b*us; 
        if SOLVER.AXI==0
            E=zeros(dims,1);
            E(1:2)=E_(1:2);
            E(4)=E_(3);
        else
            E=E_;
        end
        
        for i=1:dims
            Mat_state.Es((e-1)*dims+i,1)=E(i,1);
        end
                
        if SOLVER.UW==1 || SOLVER.UW==4
            ndw = MAT_POINT{2}(e).near;
            nn2 = length(ndw);         
            b2  = MAT_POINT{2}(e).B;

            uw=zeros(nn2*sp,1);     %Incremental displacements _ water
            for i=1:nn2
                for l=1:sp
                    uw((i-1)*sp+l,1)=d((ndw(i)-1)*df+sp+l,1);
                end
            end

            Ew_=b2*uw;
            
            if SOLVER.AXI==0
                Ew=zeros(dims,1);
                Ew(1:2)=Ew_(1:2);
                Ew(4)=Ew_(3);
            else
                Ew=Ew_;
            end

            for i=1:dims
                Mat_state.Esw((e-1)*dims+i,1)=Ew(i,1);
            end

        elseif SOLVER.UW==2
            N=MAT_POINT{2}(e).N;
            ndw = MAT_POINT{2}(e).near;
            nn2 = length(ndw);         
            b2  = MAT_POINT{2}(e).B;
            pwv=zeros(nn2,1);     %Incremental displacements _ solid
            dN=zeros(2,nn2);
            for i=1:nn2
                dN(1,i)=b2(1,(i-1)*sp+1);
                dN(2,i)=b2(2,(i-1)*sp+2);
                pwv(i,1)=d(ndw(i)*df,1);
            end
            %Pw=N'*pwv;
            Mat_state.pw(e,1)=Mat_state.pw(e,3)+N'*pwv;
            dPw=dN*pwv;
            Mat_state.dpw((e-1)*sp+1:e*sp,1)=dPw;
        end
        
        if SOLVER.REMAPPING
            %Principal strecht
            [EP_]=Principal(E_);
            MAT_POINT{1}(e).EP(:,1)=EP_(:);
        end

    end        
end

function [EP]=Principal(C)

    [eps]=eig(C);

    [EP(1),i]=max(eps);
    eps(i)=-1e32;
    [EP(2),i]=max(eps);
    eps(i)=-1e32;
    EP(3)=max(eps);


end



