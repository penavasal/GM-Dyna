
function [Mat_state,MAT_POINT]=update_F(d,Mat_state,MAT_POINT)
    
    global GEOMETRY SOLVER
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    dimf=GEOMETRY.f_dim;
    
    [aa,bb]=size(GEOMETRY.patch_con);
    area_p_n=zeros(aa,1);
    area_p_w=zeros(aa,1);
    logJ=zeros(aa,1);
    logJ_w=zeros(aa,1);

        
    m_I=AUX.m2v(eye(3));
    f_v=zeros(dimf,1);
    f_v_w=zeros(dimf,1);

    jacobians=zeros(GEOMETRY.mat_points,1);
        
    % Update **********************
    for j=1:aa
        for k=1:bb
            e=GEOMETRY.patch_con(j,k);
            nd = MAT_POINT(e).near;
            nn = length(nd);
            

            b=MAT_POINT(e).B;
            if SOLVER.AXI
                bs=zeros(5,nn*2);
                bw=zeros(5,nn*2);
                for i=2:2:nn*2
                    bs(1,i-1)=b(1,i-1);
                    bs(2,i-1)=b(2,i);
                    bs(3,i)=b(3,i);
                    bs(4,i)=b(3,i-1);
                    bs(5,i-1)=b(4,i-1);

                    bw(1,i-1)=bs(1,i-1);
                    bw(4,i)=bs(4,i);
                    bw(5,i-1)=bs(5,i-1);
                end
                
            else
                bs=zeros(4,nn*2);
                bw=zeros(4,nn*2);
                for i=2:2:nn*2
                    bs(1,i-1)=b(1,i-1);
                    bs(2,i-1)=b(2,i);
                    bs(3,i)=b(3,i);
                    bs(4,i)=b(3,i-1);

                    bw(1,i-1)=bs(1,i-1);
                    bw(4,i)=bs(4,i);
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
                dF_(dimf)=dF_(1);
            end
            dF_=dF_+m_I;

            %Vector F to Matrix F
            for i=1:5
                f_v(i,1)=Mat_state.F((e-1)*5 + i,2);
            end           
            [F]=AUX.v2m(f_v);            
            [dF]=AUX.v2m(dF_);
            F=dF*F;           
            jacobians(e)=det(F);
            
            
            if SOLVER.UW==1
                uw=zeros(nn*sp,1);     %Incremental displacements _ water
                for i=1:nn
                    for l=1:sp
                        uw((i-1)*sp+l,1)=d((nd(i)-1)*df+sp+l,1)...
                            -d((nd(i)-1)*df+sp+l,2);
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
                [F_w]=AUX.v2m(f_v_w);
                [dFw]=AUX.v2m(dF_w);
                F_w=dFw*F_w;
                
            end
            
            if SOLVER.F_BAR>0
                %%%%%%   F_BAR   %%%%%%
                F_store(k)={F};
                % New vol, dens, jaco
                area_p_n(j)=area_p_n(j)+GEOMETRY.Area(e)*jacobians(e);
                logJ(j)=logJ(j) + GEOMETRY.Area(e)*jacobians(e)*log(jacobians(e));
            else
                % New jaco
                MAT_POINT(e).J=jacobians(e);
                
                %Storage of vector F
                [f]=AUX.m2v(F);
                for i=1:dimf
                    Mat_state.F(e*dimf+1-i,1)=f(dimf+1-i);
                end

                if isnan(F)
                    disp('Fail in F');
                end
                
                if SOLVER.REMAPPING
                    %Principal strecht
                    C=F_'*F_;
                    [EP_]=Principal(C);
                    MAT_POINT(e).EP(:,1)=EP_(:);
                end
            end
            if SOLVER.UW==1
                if SOLVER.F_BAR_W>0
                    %%%%%%   F_BAR W  %%%%%%
                    F_store_w(k)={F_w};
                    % New vol, dens, jaco
                    area_p_w(j)=area_p_w(j)+GEOMETRY.Area(e)*det(F_w);
                    logJ_w(j)=logJ_w(j) + GEOMETRY.Area(e)*det(F_w)*log(det(F_w));
                else
                    %Storage of vector F
                    [f_w]=AUX.m2v(F_w);
                    for i=1:dimf
                        Mat_state.Fw(e*dimf+1-i,1)=f_w(dimf+1-i);
                    end
                end
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
                MAT_POINT(e).J=jacobians(e);
                
                %Storage of vector F
                [f]=AUX.m2v(F_);
                for i=1:dimf
                    Mat_state.F(e*dimf+1-i,1)=f(dimf+1-i);
                end
                
                if SOLVER.REMAPPING
                    %Principal strecht
                    C=F_'*F_;
                    [EP_]=Principal(C);
                    MAT_POINT(e).EP(:,1)=EP_(:);
                end
            end
        end
        if SOLVER.UW==1
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
                    [f_w]=AUX.m2v(F_w_);
                    for i=1:dimf
                        Mat_state.Fw(e*dimf+1-i,1)=f_w(dimf+1-i);
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        
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



