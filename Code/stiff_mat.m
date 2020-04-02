
function [KG]=stiff_mat(MAT_POINT,Mat_state,e,KG,T,A,BLCK)

    global SOLVER
    
    if SOLVER.SMALL==0
        [KG]=stiff_mat_L(MAT_POINT,Mat_state,e,KG,T,A,BLCK);
    else
        [KG]=stiff_mat_S(MAT_POINT,Mat_state,e,KG,A,BLCK);
    end

end

function [KG]=stiff_mat_S(MAT_POINT,Mat_state,e,KG,A,BLCK)

    global GEOMETRY SOLVER

    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
    
    % Derivatives
    nb=MAT_POINT{1}(e).near;
    b=MAT_POINT{1}(e).B;
    n =length(nb);
    
    if SOLVER.UW>0
        nbw=MAT_POINT{2}(e).near;
        bw=MAT_POINT{2}(e).B;
        nw =length(nbw);
    end
    
    volume=GEOMETRY.Area(e)*MAT_POINT{1}(e).J;
    if SOLVER.AXI
        dim=4;
        if SOLVER.UW==0 || SOLVER.UW==2
            vol=-2*pi*MAT_POINT{1}(e).xg(1)*volume; % Por que negativo??
        elseif SOLVER.UW==1
            vol=-2*pi*MAT_POINT{1}(e).xg(1)*volume;
        end
        
    else
        dim=3;
        vol=-volume;
    end
    
    % ----------------------------
    % Elemental matrices
    % ----------------------------
    D=zeros(dim);
    D(1:dim,1:dim)=A(1:dim,1:dim);
    % Material
    K_mat=b'*D*b;    
    
    if SOLVER.UW==0
        K_el=vol*K_mat;
        [KG]=assemble(KG,K_el,nb,nb,n,n,df);
    elseif SOLVER.UW==2
        Nw= MAT_POINT{2}(e).N;
        [KG]=Mat_UPw(...
            KG,vol,b,bw,nb,nbw,n,nw,Nw,K_mat,Mat_state.k(e),sp,df);      
    elseif SOLVER.UW==1
        % Material
        [KG]=Mat_UW(KG,vol,K_mat,b,bw,nb,nbw,n,nw,e,MAT_POINT{1}(e).J,sp,df,BLCK);
    end
     
end

function [KG]=stiff_mat_L(MAT_POINT,Mat_state,e,KG,T,A,BLCK)

    global GEOMETRY SOLVER

    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
    
    % Derivatives
    nb=MAT_POINT{1}(e).near;
    b=MAT_POINT{1}(e).B;
    n =length(nb);
    
    if SOLVER.UW>0
        nbw=MAT_POINT{2}(e).near;
        bw=MAT_POINT{2}(e).B;
        nw =length(nbw);
    end
    
    volume=GEOMETRY.Area(e)*MAT_POINT{1}(e).J;
    if SOLVER.AXI
        dim=4;
        if SOLVER.UW==0 || SOLVER.UW==2
            vol=-2*pi*MAT_POINT{1}(e).xg(1)*volume; % Por que negativo??
        elseif SOLVER.UW==1
            vol=-2*pi*MAT_POINT{1}(e).xg(1)*volume;
        end
        
    else
        dim=3;
        if SOLVER.UW==0 || SOLVER.UW==2
            vol=-volume;   % Por que negativo??
        elseif SOLVER.UW==1
            vol=-volume;
        end 
    end
    
    % ----------------------------
    % Elemental matrices
    % ----------------------------
    D=zeros(dim);
    D(1:dim,1:dim)=A(1:dim,1:dim);
    % Material
    K_mat=b'*D*b;    
    % Geometrical
    [K_geo]=Geo(b,T,n,sp);
    
    if SOLVER.UW==0
        K_el=vol*(K_mat+K_geo);
        [KG]=assemble(KG,K_el,nb,nb,n,n,df);
    elseif SOLVER.UW==2
        Nw= MAT_POINT{2}(e).N;
        [KG]=Mat_UPw(...
            KG,vol,b,bw,nb,nbw,n,nw,Nw,K_mat+K_geo,Mat_state.k(e),sp,df);      
    elseif SOLVER.UW==1
        % Material
        [KG]=Mat_UW(KG,vol,K_mat,b,bw,nb,nbw,n,nw,e,MAT_POINT{1}(e).J,sp,df,BLCK);
        % Geometrical
        [KG]=Geo_UW(KG,vol,e,Mat_state,MAT_POINT,b,n,nb,T,sp,df,BLCK);
    end
     
end

function [KG]=assemble(KG,K_el,nb1,nb2,n1,n2,df)
    % ----------------------------
    % Assembling of K
    % ----------------------------
    for j=1:n1
        for l=1:n2
            for k=1:df
                for r=1:df
                    KG(nb1(j)*df+1-k,nb2(l)*df+1-r)...
                    =KG(nb1(j)*df+1-k,nb2(l)*df+1-r)...
                    - K_el(j*df+1-k,l*df+1-r);
                end
            end
        end
    end
end

function [KG]=Mat_UW(KG,vol,Kel,b,bw,nb,nbw,n,nw,i,J,sp,df,BLCK)

    global SOLVER MATERIAL GEOMETRY
    
    Material=GEOMETRY.material;
    MAT=MATERIAL(BLCK).MAT;
    
    K_w=MAT{28,Material(i)};
    K_s=MAT{27,Material(i)};

    if SOLVER.AXI
        m=[1 1 0 1];
    else
        m=[1 1 0];
    end
    nn=1-(1-MAT{16,Material(i)})/J;
    Q=1/(nn/K_w+(1-nn)/K_s);

    dN=m*b;
    dNw=m*bw;
    
    Qu=Q*vol*(dN')*dN;
    Quw=Q*vol*dN'*dNw;
    Qwu=Quw';
    Qw=Q*vol*(dNw')*dNw;
    
    K_mat=vol*Kel + Qu;
    
    K_11=zeros(df*n);
    K_22=zeros(df*nw);
    K_12=zeros(df*n,df*nw);
    K_21=zeros(df*nw,df*n);
    for j=1:n
        for k=1:n
            K_11((j-1)*df+1:(j-1)*df+sp,(k-1)*df+1:(k-1)*df+sp)=...
                -K_mat((j-1)*sp+1:j*sp,(k-1)*sp+1:k*sp);
        end
    end
    [KG]=assemble(KG,K_11,nb,nb,n,n,df);
    for j=1:nw
        for k=1:nw
            K_22((j-1)*df+sp+1:j*df,(k-1)*df+sp+1:k*df)=...
                -Qw((j-1)*sp+1:j*sp,(k-1)*sp+1:k*sp);
        end
    end
    [KG]=assemble(KG,K_22,nbw,nbw,nw,nw,df);
    for j=1:n
        for k=1:nw
            K_12((j-1)*df+1:(j-1)*df+sp,(k-1)*df+sp+1:k*df)=...
                -Quw((j-1)*sp+1:j*sp,(k-1)*sp+1:k*sp);
        end
    end
    [KG]=assemble(KG,K_12,nb,nbw,n,nw,df);
    for j=1:nw
        for k=1:n
            K_21((j-1)*df+sp+1:j*df,(k-1)*df+1:(k-1)*df+sp)=...
                -Qwu((j-1)*sp+1:j*sp,(k-1)*sp+1:k*sp);
        end
    end
    [KG]=assemble(KG,K_21,nbw,nb,nw,n,df);
end

function [KG]=Mat_UPw(KG,vol,b,bw,nb,nbw,n,nw,Nw,K,k,sp,df)

    global SOLVER

    if SOLVER.AXI
        m=[1 1 0 1];
    else
        m=[1 1 0];
    end
    
    Q=b'*m'*Nw';
    
    dNw=zeros(2,nw);
    for j=1:nw
        dNw(1,j)=bw(1,(j-1)*sp+1);
        dNw(2,j)=bw(2,(j-1)*sp+2);
    end
    
    H=k*(dNw'*dNw);
  
    % ----------------------------
    % Material stiffness
    % ----------------------------
    K_1=zeros(df*n);
    K_2=zeros(df*nw);
    K_3=zeros(df*n,df*nw);
    for j=1:n
        for k=1:n
            K_1((j-1)*df+1:j*df-1,(k-1)*df+1:k*df-1)=...
                -vol*K((j-1)*sp+1:j*sp,(k-1)*sp+1:k*sp);
        end
    end
    [KG]=assemble(KG,K_1,nb,nb,n,n,df);
    for j=1:nw
        for k=1:nw
            K_2(j*df,k*df)=-vol*H(j,k);
        end
    end
    [KG]=assemble(KG,K_2,nbw,nbw,nw,nw,df);
    for j=1:n
        for k=1:nw
            K_3((j-1)*df+1:j*df-1,k*df)=...
                +vol*Q((j-1)*sp+1:j*sp,k);
        end
    end
    [KG]=assemble(KG,K_3,nb,nbw,n,nw,df);
end

function [KG]=Geo_UW(KG,vol,e,Mat_state,MAT_POINT,b,n,nb,T,sp,df,BLCK)

    global SOLVER MATERIAL GEOMETRY
    
    Material=GEOMETRY.material;
    MAT=MATERIAL(BLCK).MAT;
    
    pw=Mat_state.pw(e,1);
    % ----------------------------
    % Geometric stiffness
    % ----------------------------
    if SOLVER.AXI
        dim=sp*2+1;
        sh=zeros(n*2,5);
        for j=1:n
            sh(j*2-1,1)=b(1,j*2-1);
            sh(j*2-1,2)=b(2,j*2);
            sh(j*2,3)=b(1,j*2-1);
            sh(j*2,4)=b(2,j*2);
            sh(j*2-1,5)=b(4,j*2-1);
        end
    else
        dim=sp*2;
        sh=zeros(n*2,4);
        for j=1:n
            sh(j*2-1,1)=b(1,j*2-1);
            sh(j*2-1,2)=b(2,j*2);
            sh(j*2,3)=b(1,j*2-1);
            sh(j*2,4)=b(2,j*2);
        end
    end

    Sigma=zeros(dim);
    Sigma(1:2,1:2)=T(1:2,1:2);
    Sigma(3:4,3:4)=T(1:2,1:2);
    if SOLVER.AXI
        Sigma(5,5)=T(3,3);
    end
    
    %Water
    nn=1-(1-MAT{16,Material(e)})/MAT_POINT(e).J;
    [M4]=linearization(nn,b,nb,sp);
    
    %Solid
    K_geo_=sh*Sigma*sh';
    K_geo(n*df,n*df)=0;
    for j=1:n
        for r=1:n
            for k=1:sp
                for l=1:sp
                    K_geo((j-1)*df+k,(r-1)*df+l)=K_geo((j-1)*df+k,(r-1)*df+l)+...
                        K_geo_((j-1)*sp+k,(r-1)*sp+l);
                    K_geo((j-1)*df+k,(r-1)*df+l)=K_geo((j-1)*df+k,(r-1)*df+l)+...
                        pw*M4((j-1)*sp+k,(r-1)*sp+l);
                    K_geo((j-1)*df+k+2,(r-1)*df+l)=K_geo((j-1)*df+k+2,(r-1)*df+l)+...
                        pw*M4((j-1)*sp+k,(r-1)*sp+l);
                end
            end
        end
    end
    [KG]=assemble(KG,K_geo*vol,nb,nb,n,n,df);
end

function [K_geo]=Geo(b,T,n,sp)

    global SOLVER
    % ----------------------------
    % Geometric stiffness
    % ----------------------------
    if SOLVER.AXI
        dim=sp*2+1;
        sh=zeros(n*2,5);
        for j=1:n
            sh(j*2-1,1)=b(1,j*2-1);
            sh(j*2-1,2)=b(2,j*2);
            sh(j*2,3)=b(1,j*2-1);
            sh(j*2,4)=b(2,j*2);
            sh(j*2-1,5)=b(4,j*2-1);
        end
    else
        dim=sp*2;
        sh=zeros(n*2,4);
        for j=1:n
            sh(j*2-1,1)=b(1,j*2-1);
            sh(j*2-1,2)=b(2,j*2);
            sh(j*2,3)=b(1,j*2-1);
            sh(j*2,4)=b(2,j*2);
        end
    end

    Sigma=zeros(dim);
    Sigma(1:2,1:2)=T(1:2,1:2);
    Sigma(3:4,3:4)=T(1:2,1:2);
    if SOLVER.AXI
        Sigma(5,5)=T(3,3);
    end
    
    K_geo=sh*Sigma*sh';
end

function [M4]=linearization(n,b,nb,sp)

    global SOLVER

    m=length(nb);
    sh=zeros(m*sp,sp);
    v=zeros(m*sp,sp);
    N=zeros(m*sp,1);
    for j=1:m
        sh((j-1)*sp+1,1)=b(1,(j-1)*sp+1);
        sh((j-1)*sp+1,2)=b(2,(j-1)*sp+2);
        
        N((j-1)*sp+1,1)=b(1,(j-1)*sp+1);
        N((j-1)*sp+2,1)=b(2,(j-1)*sp+2);
        
        if SOLVER.AXI
            sh((j-1)*sp+1,1)=sh((j-1)*sp+1,1)+b(4,(j-1)*sp+1);
            N((j-1)*sp+1,1) = N((j-1)*sp+1,1)+b(4,(j-1)*sp+1);
        end
        
        sh((j-1)*sp+2,1)=sh((j-1)*sp+1,1);
        sh((j-1)*sp+2,2)=sh((j-1)*sp+1,2);
        
        v((j-1)*sp+1,1)=1;
        v((j-1)*sp+2,2)=1;
    end
    
    M1=sh*v';
    M2=M1.*M1';
    M3=N*N';
    
    M4=M3+(1-n)/n*M2;

end
