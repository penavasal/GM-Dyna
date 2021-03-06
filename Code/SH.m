
classdef SH
    
    % Shape function methods
    
    methods(Static)
        function MAT_POINT=calculation(...
            INITIAL,MAT_POINT,Disp_field,NODE_LIST)

            global GEOMETRY SOLVER LME_param

            REM_T=0;

            % INITIAL CALCULATION

            if INITIAL==1

                % Measurement of remapping
                [phases,~]=size(SOLVER.PHASES);
                for k=1:phases
                    for i=1:GEOMETRY.mat_points
                        for j=1:3
                            MAT_POINT{k}(i).EP(j,1)=1;
                            MAT_POINT{k}(i).EP(j,2)=1;
                        end  
                    end  
                end

                % Only for LME shape function
                if strcmp(SOLVER.TYPE{2},'LME')

                    MAT_POINT=LME.initialize(MAT_POINT,NODE_LIST);

                    [phases,~]=size(SOLVER.PHASES);
                    [shf,~]=size(LME_param);
                    for sh=1:shf
                        for ph=1:phases
                            if SOLVER.PHASES{ph,2}==sh
                                break;
                            end
                        end

                        n_sp=LME.nodalspacing(MAT_POINT{ph});
                        for i=1:GEOMETRY.mat_points
                            xg=GEOMETRY.xg_0(i,:);
                            [MAT_POINT{ph}]=LME.near(i,GEOMETRY.x_0,xg,...
                                n_sp,MAT_POINT{ph},sh,ph); 
                            [MAT_POINT{ph},FAIL]=LME.shapef...
                                (i,GEOMETRY.x_0,n_sp,xg,MAT_POINT{ph},sh);
                            if FAIL
                                error('FAIL in the initial calculation of Shape function \n');
                            end
                        end

                        for ph1=ph+1:phases
                            if SOLVER.PHASES{ph1,2}==sh
                                for i=1:GEOMETRY.mat_points
                                    MAT_POINT{ph1}(i).N = MAT_POINT{ph}(i).N;
                                    MAT_POINT{ph1}(i).B = MAT_POINT{ph}(i).B;
                                    MAT_POINT{ph1}(i).xi= MAT_POINT{ph}(i).xi;
                                    MAT_POINT{ph1}(i).near= MAT_POINT{ph}(i).near;
                                end
                            end
                        end
                    end

                else

                    [phases,~]=size(SOLVER.PHASES);

                    if strcmp(GEOMETRY.ELEMENT,'Q8P4-4') || strcmp(GEOMETRY.ELEMENT,'Q8P4')...
                            || strcmp(GEOMETRY.ELEMENT,'T6P3-3') || strcmp(GEOMETRY.ELEMENT,'T6P3')
                        shf=2;
                        if SOLVER.UW==3
                            SOLVER.PHASES{3,2}=2;
                        else
                            SOLVER.PHASES{2,2}=2;
                        end
                    else
                        shf=1;
                    end

                    for sh=1:shf
                        for ph=1:phases
                            if SOLVER.PHASES{ph,2}==sh
                                break;
                            end
                        end

                        for i=1:GEOMETRY.mat_points

                            %e=MAT_POINT{ph}(i).element;

                            %elem=GEOMETRY.elem(e,:);
                            elem=MAT_POINT{ph}(i).near;

                            % Discover isoparametric coordinates
                            coord=zeros(length(elem),GEOMETRY.sp);
                            for j=1:length(elem)
                                coord(j,:)=GEOMETRY.x_0(elem(j),:);
                            end
                            [xilist_0,coord_xi] = ...
                                FEM.integrationpoints(GEOMETRY.sp,length(elem),1);
                            [xilist]=FEM.search_xilist_NR(coord,GEOMETRY.xg_0(i,:),...
                                xilist_0,coord_xi);

                            MAT_POINT{ph}(i).xi=xilist;

                            % Calculate initial shape function
                            [MAT_POINT{ph}]=FEM.shapef(i,coord,xilist,MAT_POINT{ph});

                        end
                        for ph1=ph+1:phases
                            if SOLVER.PHASES{ph1,2}==sh
                                for i=1:GEOMETRY.mat_points
                                    MAT_POINT{ph1}(i).N = MAT_POINT{ph}(i).N;
                                    MAT_POINT{ph1}(i).B = MAT_POINT{ph}(i).B;
                                    MAT_POINT{ph1}(i).xi= MAT_POINT{ph}(i).xi;
                                    MAT_POINT{ph1}(i).near= MAT_POINT{ph}(i).near;
                                end
                            end
                        end
                    end

                end

            else

                if strcmp(SOLVER.TYPE{2},'LME')

                    [phases,~]=size(SOLVER.PHASES);
                    [shf,~]=size(LME_param);
                    for sh=1:shf
                        for ph=1:phases
                            if SOLVER.PHASES{ph,2}==sh
                                break;
                            end
                        end

                        for i=1:GEOMETRY.mat_points            
                            [MAT_POINT{ph},REMAP]=LME.REMAPPING(MAT_POINT{ph},i);
                            REM_T=max(REM_T,REMAP);  
                            if REMAP
                                if SOLVER.TYPE{1}==0 %OTM
                                    x=Disp_field.x_a;
                                    jacobians=MAT_POINT{1}(i).J;
                                    h=GEOMETRY.h_ini.*sqrt(jacobians);
                                    xg = getfield(MAT_POINT(i),'xg');
                                elseif SOLVER.TYPE{1}==1 %MPM
                                    x=GEOMETRY.x_0;
                                    xg=GEOMETRY.xg_0(i,:);
                                    h=GEOMETRY.h_ini;
                                end

                                [MAT_POINT{ph}]=LME.near(i,x,xg,h,MAT_POINT{ph}); 
                                [MAT_POINT{ph},FAIL]=LME.shapef(i,x,h,xg,MAT_POINT{ph});
                                if FAIL
                                    error('FAIL in the initial calculation of Shape function \n');
                                end
                            end
                        end
                        for ph1=ph+1:phases
                            if SOLVER.PHASES{ph1,2}==sh
                                for i=1:GEOMETRY.mat_points
                                    MAT_POINT{ph1}(i).N = MAT_POINT{ph}(i).N;
                                    MAT_POINT{ph1}(i).B = MAT_POINT{ph}(i).B;
                                    MAT_POINT{ph1}(i).xi= MAT_POINT{ph}(i).xi;
                                    MAT_POINT{ph1}(i).near= MAT_POINT{ph}(i).near;
                                end
                            end
                        end
                    end
                else

                    if SOLVER.TYPE{1}==1 %MPM

                        for i=1:GEOMETRY.mat_points

                            el=MAT_POINT(i).element;
                            [I]=IoO(MAT_POINT(i).xg,...
                                Disp_field.x_a,GEOMETRY.elem(el,:));


                            if I==0 || el==0
                                el_near=GEOMETRY.element_near{e};
                                j=0;
                                while I==0
                                    j=j+1;
                                    e=el_near(j);
                                    if j>length(el_near)
                                        disp('Please, reduce tim step, mp_s jump 2 grid elements')
                                        stop
                                    end
                                    if e~=el
                                        [I]=IoO(MAT_POINT(i).xg,Disp_field.x_a,...
                                            GEOMETRY.elem(e,:));  
                                    end
                                end


                                elem=GEOMETRY.elem(e,:);
                                % Discover isoparametric coordinates
                                coord=zeros(length(elem),GEOMETRY.sp);
                                for j=1:length(elem)
                                    coord(j,:)=GEOMETRY.x_0(elem(j),:);
                                end
                                [xilist_0,coord_xi] = ...
                                    FEM.integrationpoints(GEOMETRY.sp,length(elem),1);
                                [xilist]=FEM.search_xilist_NR(coord,MAT_POINT(i).xg,...
                                    xilist_0,coord_xi);
                            else
                                elem=GEOMETRY.elem(el,:);
                                % Discover isoparametric coordinates
                                coord=zeros(length(elem),GEOMETRY.sp);
                                for j=1:length(elem)
                                    coord(j,:)=GEOMETRY.x_0(elem(j),:);
                                end
                                [~,coord_xi] = ...
                                    FEM.integrationpoints(GEOMETRY.sp,length(elem),1);
                                [xilist]=FEM.search_xilist_NR(coord,MAT_POINT(i).xg,...
                                    MAT_POINT(i).xi,coord_xi);
                                %[xilist]=FEM.search_xilist(coord,MAT_POINT(i).xg,xilist,i);
                            end

                            MAT_POINT(i).xi=xilist;

                            % Calculate initial shape function
                            [MAT_POINT]=FEM.shapef(i,coord,xilist,MAT_POINT);
                        end

                    else %FEM
                        REM_T=1;
                        [phases,~]=size(SOLVER.PHASES);

                        if strcmp(GEOMETRY.ELEMENT,'Q8P4-4') || strcmp(GEOMETRY.ELEMENT,'Q8P4')...
                                || strcmp(GEOMETRY.ELEMENT,'T6P3-3') || strcmp(GEOMETRY.ELEMENT,'T6P3')
                            shf=2;
                            if SOLVER.UW==3
                                SOLVER.PHASES{3,2}=2;
                            else
                                SOLVER.PHASES{2,2}=2;
                            end
                        else
                            shf=1;
                        end

                        for sh=1:shf
                            for ph=1:phases
                                if SOLVER.PHASES{ph,2}==sh
                                    break;
                                end
                            end

                            for i=1:GEOMETRY.mat_points
                                
                                elem=MAT_POINT{ph}(i).near;
                                dis=MAT_POINT{ph}(i).REM;
                                Dis=max(abs(dis));
                                
                                if Dis>0

                                    % Discover isoparametric coordinates
                                    coord=zeros(length(elem),GEOMETRY.sp);
                                    for j=1:length(elem)
                                        coord(j,:)=Disp_field.x_a(elem(j),:);
                                    end

                                    [~,coord_xi] = ...
                                        FEM.integrationpoints(GEOMETRY.sp,length(elem),1);
                                    [xilist]=FEM.search_xilist_NR(coord,MAT_POINT{ph}(i).xg,...
                                        MAT_POINT{ph}(i).xi,coord_xi);

                                    MAT_POINT{ph}(i).xi=xilist;


                                    % Calculate initial shape function
                                    [MAT_POINT{ph}]=FEM.shapef(i,coord,xilist,MAT_POINT{ph});
                                end
                            end
                            for ph1=ph+1:phases
                                if SOLVER.PHASES{ph1,2}==sh
                                    for i=1:GEOMETRY.mat_points
                                        MAT_POINT{ph1}(i).N = MAT_POINT{ph}(i).N;
                                        MAT_POINT{ph1}(i).B = MAT_POINT{ph}(i).B;
                                        MAT_POINT{ph1}(i).xi= MAT_POINT{ph}(i).xi;
                                        MAT_POINT{ph1}(i).near= MAT_POINT{ph}(i).near;
                                    end
                                end
                            end
                        end
                    end
                end
            end

            if REM_T || INITIAL
                if SOLVER.TYPE{1}==0 %OTM
                    [jacobians]=LIB.S2list(MAT_POINT{1},'J');
                    volume=GEOMETRY.Area.*jacobians;
                elseif SOLVER.TYPE{1}==1 || INITIAL %MPM || 
                    volume=GEOMETRY.Area;
                end

                % MAKE Bbar - Patch
                for ph=1:phases
                    if SOLVER.B_BAR==1
                        [MAT_POINT{ph}]=SH.bbar_matrix(GEOMETRY.patch_con,GEOMETRY.patch_el,...
                            GEOMETRY.mat_points,MAT_POINT{ph},volume,SOLVER.B_BAR);
                    else
                        [MAT_POINT{ph}]=SH.b_m(GEOMETRY.mat_points,GEOMETRY.sp,MAT_POINT{ph},INITIAL);
                    end
                end
            end

        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % B and B-bar matrix construction functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [MAT_POINT]=b_m(elements,sp,MAT_POINT,INIT)

            global SOLVER

            if SOLVER.AXI
                for i=1:elements
                    r=MAT_POINT(i).xg(1);
                    nb=MAT_POINT(i).near;
                    sh=MAT_POINT(i).B;
                    pp=MAT_POINT(i).N;
                    n=length(nb);
                    b2=zeros(4,sp*n);
                    for j=1:n
                        b2(1,sp*j-1)=sh(j,1);
                        b2(2,sp*j)  =sh(j,2);
                        b2(4,sp*j-1)=pp(j)/r;
                        b2(3,sp*j-1)=sh(j,2);
                        b2(3,sp*j)  =sh(j,1);
                    end
                    MAT_POINT(i).B=b2;
                    clear b2;
                end
            else
                for i=1:elements 
                    dis=MAT_POINT(i).REM;
                    Dis=max(abs(dis));
                                
                    if Dis>0 || INIT==1
                        nb=MAT_POINT(i).near;
                        sh=MAT_POINT(i).B;
                        n=length(nb);
                        b=zeros(3,sp*n);
                        for j=1:n
                            b(1,sp*j-1)=sh(j,1);
                            b(2,sp*j)  =sh(j,2);
                            b(3,sp*j-1)=sh(j,2);
                            b(3,sp*j)  =sh(j,1);
                        end
                        MAT_POINT(i).B=b;
                        clear b;
                    end
                end
            end

        end

        %function [BBar2,near,p]=bbar_matrix(x_a,patch_con,patch_el,elem,near,...
        %   dp,p,Area,BBAR)

        function[MAT_POINT]=bbar_matrix...
                        (patch_con,patch_el,elements,MAT_POINT,Area,BBAR)

            global SOLVER

            [MAT_POINT]=SH.new_nb(MAT_POINT,patch_con,elements);
            if SOLVER.AXI==0
                [MAT_POINT]=SH.bbm(patch_con,patch_el,MAT_POINT,Area,BBAR);
            elseif SOLVER.AXI==0
                [MAT_POINT]=SH.bbm_axi(patch_con,patch_el,MAT_POINT,Area,BBAR);
            end

        end

        function [MAT_POINT]=bbm(patch_con,patch_el,MAT_POINT,Area,BBAR)

            global GEOMETRY SOLVER

            [patch,EN]=size(patch_con);
            [elements,~]=size(patch_el);
            df=GEOMETRY.df;
            sp=GEOMETRY.sp;

            for i=1:elements
                nb = MAT_POINT(i).near;
                sh = MAT_POINT(i).B;
                n=length(nb);
                if SOLVER.UW==0 || SOLVER.UW==2
                    b=zeros(3,sp*n);
                    t=zeros(1,n*df);
                    for j=1:n
                        for k=1:sp
                            t(j*sp+1-k)=sh(j,sp+1-k);
                        end
                        b(1,sp*j-1)=sh(j,1);
                        b(2,sp*j)  =sh(j,2);
                        b(3,sp*j-1)=sh(j,2);
                        b(3,sp*j)  =sh(j,1);
                    end
                    T(i)={t};
                    B(i)={b};
                    clear t b
                elseif SOLVER.UW==1 || SOLVER.UW==4
                    b2=zeros(5,df*n);
                    ts=zeros(1,n*df);
                    tw=zeros(1,n*df);
                    for j=1:n
                        for k=1:sp
                            ts(j*df-1-k)=sh(j,sp+1-k);
                            tw(j*df+1-k)=sh(j,sp+1-k);
                        end
                        b2(1,df*j-3)=sh(j,1);
                        b2(2,df*j-2)=sh(j,2);
                        b2(3,df*j-3)=sh(j,2);
                        b2(3,df*j-2)=sh(j,1);
                        b2(4,df*j-1)=sh(j,1);
                        b2(5,df*j)=sh(j,2);
                    end
                    Tw(i)={tw};
                    Ts(i)={ts};
                    B2(i)={b2};
                    clear b2 ts tw;
                end  
            end

            if SOLVER.UW==0 || SOLVER.UW==2
                for i=1:patch
                    t=T{patch_con(i,1)};
                    l=length(t);
                    ts(l)=0;
                    A=Area(patch_con(i,1));
                    ts(:)=A*t(:);        
                    for j=2:EN
                        t=T{patch_con(i,j)};
                        a=Area(patch_con(i,j));
                        ts(:)=ts(:)+a*t(:);
                        A=A+a;
                    end
                    Tp(i)={ts};
                    At(i)=A;
                    clear ts;
                end

                m=[1 1 0];
                for i=1:elements
                    b =B{i};
                    t =T{i};
                    A =At(patch_el(i));
                    tt=Tp{patch_el(i)};
                    bb=b-1/sp*m'*(t-1/A*tt);
                    MAT_POINT(i).B=bb;
                    clear b t tt bb;
                end
            elseif SOLVER.UW==1 || SOLVER.UW==4
                for i=1:patch
                    tw=Tw{patch_con(i,1)};
                    ts=Ts{patch_con(i,1)};
                    l=length(ts);
                    l2=length(tw);
                    tls(l)=0;
                    tlw(l2)=0;
                    A=Area(patch_con(i,1));
                    tls(:)=A*ts(:);
                    tlw(:)=A*tw(:);
                    for j=2:EN
                        ts=Ts{patch_con(i,j)};
                        tw=Tw{patch_con(i,j)};
                        a=Area(patch_con(i,j));
                        tls(:)=tls(:)+a*ts(:);
                        tlw(:)=tlw(:)+a*tw(:);
                    end
                    Tps(i)={tls};
                    Tpw(i)={tlw};
                    clear tls tlw;
                end


                for i=1:elements

                if BBAR(i)
                    ms=[1 1 0 0 0];
                    mw=[0 0 0 1 1];
                else
                    ms=[0 0 0 0 0];
                    mw=[0 0 0 0 0];
                end

                    b  =B2{i};
                    tw =Tw{i};
                    ts =Ts{i};
                    A =At(patch_el(i));
                    ttw=Tpw{patch_el(i)};
                    tts=Tps{patch_el(i)};
                    bb=b-1/sp*ms'*(ts-1/A*tts)-1/sp*mw'*(tw-1/A*ttw);
                    MAT_POINT(i).B=bb;
                    clear b tw ts tts ttw bb;
                end
            end
        end

        function [MAT_POINT]=new_nb(MAT_POINT,patch_con,elements)

            [patch,EN]=size(patch_con);

            for i=1:patch
                nw_near_=MAT_POINT(patch_con(i,1)).near;
                k=length(nw_near_);
                for j=2:EN
                    near_=MAT_POINT(patch_con(i,j)).near;
                    s=length(near_);
                    for l=1:s
                        nd=near_(l);
                        t=0;
                        for r=1:k
                            if nd==nw_near_(r)
                                t=1;
                            end
                        end
                        if t==0
                            k=k+1;
                            nw_near_(k)=nd;
                        end
                    end
                end
                for j=1:EN
                    el=patch_con(i,j);
                    nw_near(el)={nw_near_};
                end
            end

            for i=1:elements
                nw_near_= nw_near{i};
                near_   = MAT_POINT(i).near;
                p_      = MAT_POINT(i).N;
                dp_     = MAT_POINT(i).B;

                m =length(nw_near_);
                mm=length(near_);

                for j=1:m
                    nd=nw_near_(j);
                    t=0;
                    for k=1:mm
                        if nd==near_(k)
                            t=k;
                        end
                    end
                    if t
                        nw_p_(j,1) = p_(t,1);
                        nw_dp_(j,1)= dp_(t,1);
                        nw_dp_(j,2)= dp_(t,2);
                    else
                        nw_p_(j,1) = 0;
                        nw_dp_(j,1)= 0;
                        nw_dp_(j,2)= 0;
                    end
                end
               MAT_POINT(i).N=nw_p_; 
               MAT_POINT(i).B=nw_dp_;
               MAT_POINT(i).near=nw_near_;
            end
        end

        function [MAT_POINT]=bbm_axi(patch_con,patch_el,MAT_POINT,Area,BBAR)

            global GEOMETRY SOLVER

            [patch,EN]=size(patch_con);
            [elements,~]=size(patch_el);

            df=GEOMETRY.df;
            sp=GEOMETRY.sp;
            bdim=GEOMETRY.b_dim;

            sp1=sp+1;
            if SOLVER.UW==0 || SOLVER.UW==2
                for i=1:elements
                    r=MAT_POINT(i).xg(1);
                    nb = MAT_POINT(i).near;
                    sh = MAT_POINT(i).B;
                    pp = MAT_POINT(i).N;
                    n=length(nb);
                    b=zeros(bdim,sp*n);
                    t=zeros(1,n*sp);
                    for j=1:n
                        for k=1:sp
                            t(j*sp+1-k)=sh(j,sp+1-k);
                        end
                        t(j*sp-1)=t(j*sp-1)+pp(j)/r;
                        b(1,sp*j-1)=sh(j,1);
                        b(2,sp*j)  =sh(j,2);
                        b(3,sp*j-1)=pp(j)/r;
                        b(4,sp*j-1)=sh(j,2);
                        b(4,sp*j)  =sh(j,1);
                    end
                    T(i)={t};
                    B(i)={b};
                    clear t b;
                end

                for i=1:patch
                    t=T{patch_con(i,1)};
                    l=length(t);
                    ts(l)=0;
                    A=Area(patch_con(i,1));
                    ts(:)=A*t(:);        
                    for j=2:EN
                        t=T{patch_con(i,j)};
                        a=Area(patch_con(i,j));
                        ts(:)=ts(:)+a*t(:);
                        A=A+a;
                    end
                    Tp(i)={ts};
                    At(i)=A;
                    clear ts;
                end

                if BBAR
                    m=[1 1 1 0];
                else
                    m=[0 0 0 0];
                end

                for i=1:elements
                    b =B{i};
                    t =T{i};
                    A =At(patch_el(i));
                    tt=Tp{patch_el(i)};
                    bb=b-1/sp1*m'*(t-1/A*tt);
                    MAT_POINT(i).B=bb;
                    clear b t tt bb;
                end
            elseif SOLVER.UW==1 || SOLVER.UW==4
                for i=1:elements
                    r=MAT_POINT(i).xg(1);
                    nb = MAT_POINT(i).near;
                    sh = MAT_POINT(i).B;
                    pp = MAT_POINT(i).N;
                    n=length(nb);
                    b2=zeros(7,df*n);
                    ts=zeros(1,n*df);
                    tw=zeros(1,n*df);
                    for j=1:n
                        for k=1:sp
                            ts(j*df-1-k)=sh(j,sp+1-k);
                            tw(j*df+1-k)=sh(j,sp+1-k);
                        end
                        ts(j*df-3)=ts(j*df-3)+pp(j)/r;
                        tw(j*df-1)=tw(j*df-1)+pp(j)/r;
                        b2(1,df*j-3)=sh(j,1);
                        b2(2,df*j-2)=sh(j,2);
                        b2(3,df*j-3)=pp(j)/r;
                        b2(4,df*j-3)=sh(j,2);
                        b2(4,df*j-2)=sh(j,1);
                        b2(5,df*j-1)=sh(j,1);
                        b2(6,df*j)=sh(j,2);
                        b2(7,df*j-1)=pp(j)/r;
                    end
                    Tw(i)={tw};
                    Ts(i)={ts};
                    B2(i)={b2};
                    clear b2 ts tw;
                end

                for i=1:patch
                    A=Area(patch_con(i,1));      
                    for j=2:EN
                        a=Area(patch_con(i,j));
                        A=A+a;
                    end
                    At(i)=A;
                end

                for i=1:patch
                    tw=Tw{patch_con(i,1)};
                    ts=Ts{patch_con(i,1)};
                    l=length(ts);
                    l2=length(tw);
                    tls(l)=0;
                    tlw(l2)=0;
                    A=Area(patch_con(i,1));
                    tls(:)=A*ts(:);
                    tlw(:)=A*tw(:);
                    for j=2:EN
                        ts=Ts{patch_con(i,j)};
                        tw=Tw{patch_con(i,j)};
                        a=Area(patch_con(i,j));
                        tls(:)=tls(:)+a*ts(:);
                        tlw(:)=tlw(:)+a*tw(:);
                    end
                    Tps(i)={tls};
                    Tpw(i)={tlw};
                    clear tls tlw;
                end

               if BBAR
                    ms=[1 1 1 0 0 0 0];
                    mw=[0 0 0 0 1 1 1];
               else
                    ms=[0 0 0 0 0 0 0];
                    mw=[0 0 0 0 0 0 0];
               end


                for i=1:elements
                    b  =B2{i};
                    tw =Tw{i};
                    ts =Ts{i};
                    A =At(patch_el(i));
                    ttw=Tpw{patch_el(i)};
                    tts=Tps{patch_el(i)};
                    bb=b-1/sp1*ms'*(ts-1/A*tts)-1/sp1*mw'*(tw-1/A*ttw);
                    MAT_POINT(i).B=bb;
                    clear b tw ts tts ttw bb;
                end
            end
        end

        function MAT_POINT=remap(MAT_POINT,Disp_field)

            global SOLVER GEOMETRY
            els=GEOMETRY.mat_points;
            tol=SOLVER.REMAPPING_tol;

            REMAP=0;
            for i=1:els
                dis(3)=0;
                Ep=MAT_POINT{1}(i).EP;
                num=min(Ep(:,1),Ep(:,2));
                for j=1:3
                    dis(j)=(Ep(j,1)-Ep(j,2))/num(j);
                end
                Dis=max(abs(dis));

                if Dis>tol
                    MAT_POINT{1}(i).REM=dis;
                    REMAP=REMAP+1;
                    MAT_POINT{1}(i).EP(:,2)=Ep(:,1);
                else
                    MAT_POINT{1}(i).REM=zeros(3,1);
                end
            end
            
            if REMAP>0
                [MAT_POINT]=SH.calculation(0,MAT_POINT,Disp_field,0);
            end

        end

        % function [patch_con,patch_el]=patch_nodes(x_a,k_bd,elem)
        % 
        %     [nodes,~]=size(x_a);
        %     [bd,~]=size(k_bd);
        %     [elements,NNE]=size(elem);
        % 
        %     EN=4;
        % 
        %     nod_con=zeros(nodes,1);
        %     patch_el=zeros(elements,1);
        %     patch_nd(1)=0;
        %     patch_con(1,4)=0;
        %     p=0;
        %     for i=1:nodes
        %         l=0;
        %         for j=1:elements
        %             for k=1:NNE
        %                 if elem(j,k)==i
        %                     l=l+1;
        %                     nod_con(i,l)=j;
        %                 end
        %             end    
        %         end 
        %         b=0;
        %         for j=1:bd
        %             if k_bd(j)==i
        %                 b=1;
        %             end
        %         end
        % 
        %         %Patch Node
        %         if l==EN && b==0
        %             p=p+1;
        %             patch_nd(p)=i;
        %             for j=1:EN
        %                 patch_con(p,j)=nod_con(i,j);
        %             end
        %         end 
        %     end
        % 
        % 
        %     for j=1:p
        %         for k=1:EN
        %             patch_el(patch_con(j,k))=j;
        %         end
        %     end
        % end
    end
end
