
function [stiff_mtx]=stiff_mat(MAT_POINT,Mat_state,e,stiff_mtx,T,A)

    global GEOMETRY SOLVER

    df=GEOMETRY.df;
    sp=GEOMETRY.sp;
    
    % Derivatives
    nb=MAT_POINT(e).near;
    b=MAT_POINT(e).B;
    n =length(nb);
    
                

    
    volume=GEOMETRY.Area(e)*MAT_POINT(e).J;
    if SOLVER.AXI
        dim=4;
        if SOLVER.UW==0
            d2=dim;
            vol=-2*pi*MAT_POINT(e).xg(1)*volume; % Por que negativo??
        else
            d2=7;
            vol=2*pi*MAT_POINT(e).xg(1)*volume;
        end
        
    else
        dim=3;
        if SOLVER.UW==0
            d2=dim;
            vol=-volume;   % Por que negativo??
        else
            d2=5;
            vol=volume;
        end 
    end
    
    % ----------------------------
    % Elemental matrices
    % ----------------------------
    D=zeros(d2);
    D(1:dim,1:dim)=A(1:dim,1:dim);
    if SOLVER.UW==0
        % Material
        K_mat=b'*D*b;    
        % Geometrical
        [K_geo]=Geo(b,T,n,sp);
    elseif SOLVER.UW==1
        % Material
        [K_mat]=Mat_UW(b,n,e,D,MAT_POINT(e).J,sp,df);
        % Geometrical
        [K_geo]=Geo_UW(e,Mat_state,MAT_POINT,b,n,nb,T,sp,df);
    end
    
    
    % ----------------------------
    % Assembling of K
    % ----------------------------
    K_el=vol*(K_mat+K_geo);
    for j=1:n
        for l=1:n
            for k=1:df
                for r=1:df
                    stiff_mtx(nb(j)*df+1-k,nb(l)*df+1-r)...
                    =stiff_mtx(nb(j)*df+1-k,nb(l)*df+1-r)...
                    -K_el(j*df+1-k,l*df+1-r);
                end
            end
        end
    end


end

function [K_mat]=Mat_UW(b,n,i,D,J,sp,df)

    global SOLVER MATERIAL
    
    Material=MATERIAL.e;
    MAT=MATERIAL.MAT;
    
    K_w=MAT(28,Material(i));
    K_s=MAT(27,Material(i));

    if SOLVER.AXI
        m=[1 1 0 1 1 1 1];
    else
        m=[1 1 0 1 1];
    end
    nn=1-(1-MAT(16,Material(i)))/J;
    Q=1/(nn/K_w+(1-nn)/K_s);

    if SOLVER.AXI
        for j=1:n
            b2(1:4,(j-1)*df+1:(j-1)*df+2)=b(1:4,(j-1)*sp+1:(j-1)*sp+2);
            b2(5:6,(j-1)*df+3:(j-1)*df+4)=b(1:2,(j-1)*sp+1:(j-1)*sp+2);
            b2(7,(j-1)*df+3:(j-1)*df+4)=b(4,(j-1)*sp+1:(j-1)*sp+2);
        end
    else
        for j=1:n
            b2(1:3,(j-1)*df+1:(j-1)*df+2)=b(1:3,(j-1)*sp+1:(j-1)*sp+2);
            b2(4:5,(j-1)*df+3:(j-1)*df+4)=b(1:2,(j-1)*sp+1:(j-1)*sp+2);
        end
    end
    
    % ----------------------------
    % Material stiffness
    % ----------------------------
    D=D+Q*(m'*m);
    K_mat=b2'*D*b2;

end

function [K_geo]=Geo_UW(e,Mat_state,MAT_POINT,b,n,nb,T,sp,df)

    global SOLVER MATERIAL
    
    Material=MATERIAL.e;
    MAT=MATERIAL.MAT;
    
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
    nn=1-(1-MAT(16,Material(e)))/MAT_POINT(e).J;
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
