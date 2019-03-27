function [BBar2,near,p]=bbar_matrix(x_a,patch_con,patch_el,elem,near,...
    dp,p,Area,BBAR)

[elements,~]=size(elem);
[~,sp]=size(x_a);

[near,p,dp]=new_nb(p,dp,near,patch_con,elements);
[B,Bbar,BBar2]=bbm(patch_con,patch_el,near,dp,sp,Area,BBAR);

end


function [B,Bbar,BBar2]=bbm(patch_con,patch_el,near,dp,sp,Area,BBAR)
    [patch,EN]=size(patch_con);
    [elements,~]=size(patch_el);
    df=2*sp;

    for i=1:elements
        nb=near{i};
        sh=dp{i};
        n=length(nb);
        b=zeros(3,sp*n);
        b2=zeros(5,df*n);
        for j=1:n
            for k=1:sp
                t(j*sp+1-k)=sh(j,sp+1-k);
            end
            b(1,sp*j-1)=sh(j,1);
            b(2,sp*j)  =sh(j,2);
            b(3,sp*j-1)=sh(j,2);
            b(3,sp*j)  =sh(j,1);
        end
        ts=zeros(1,n*df);
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
        T(i)={t};
        B(i)={b};
        Tw(i)={tw};
        Ts(i)={ts};
        B2(i)={b2};
        clear t b b2 ts tw;
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

    m=[1 1 0];
    for i=1:elements
        b =B{i};
        t =T{i};
        A =At(patch_el(i));
        tt=Tp{patch_el(i)};
        bb=b-1/sp*m'*(t-1/A*tt);
        Bbar(i)={bb};
        clear b t tt bb;
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
    
    %ms=[1 1 0 1 1];
    %mw=ms;
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
        BBar2(i)={bb};
        clear b tw ts tts ttw bb;
    end
end
    
function [patch_con,patch_el]=patch_nodes(x_a,k_bd,elem)

    [nodes,~]=size(x_a);
    [bd,~]=size(k_bd);
    [elements,NNE]=size(elem);

    EN=4;

    nod_con=zeros(nodes,1);
    patch_el=zeros(elements,1);
    patch_nd(1)=0;
    patch_con(1,4)=0;
    p=0;
    for i=1:nodes
        l=0;
        for j=1:elements
            for k=1:NNE
                if elem(j,k)==i
                    l=l+1;
                    nod_con(i,l)=j;
                end
            end    
        end 
        b=0;
        for j=1:bd
            if k_bd(j)==i
                b=1;
            end
        end

        %Patch Node
        if l==EN && b==0
            p=p+1;
            patch_nd(p)=i;
            for j=1:EN
                patch_con(p,j)=nod_con(i,j);
            end
        end 
    end


    for j=1:p
        for k=1:EN
            patch_el(patch_con(j,k))=j;
        end
    end
end

function [nw_near,nw_p,nw_dp]=new_nb(p,dp,near,patch_con,elements)

    [patch,EN]=size(patch_con);

    for i=1:patch
        nw_near_=near{patch_con(i,1)};
        k=length(nw_near_);
        for j=2:EN
            near_=near{patch_con(i,j)};
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
        near_   = near{i};
        p_      = p{i};
        dp_     = dp{i};

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
       nw_p(i)={nw_p_}; 
       nw_dp(i)={nw_dp_};    
    end
end