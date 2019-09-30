function [Mat_state]=internal_forces(MAT_POINT,Mat_state)

    global GEOMETRY SOLVER
    
    sig=zeros(4,1);
    Mat_state.fint(:,1) = zeros(GEOMETRY.nodes*GEOMETRY.df,1);
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;

    for e=1:GEOMETRY.mat_points
        
        nd=MAT_POINT(e).near;
        B_=MAT_POINT(e).B;
        nn =length(nd);

        % Derivatives
        if SOLVER.AXI
            sh=zeros(3,nn);
            for i=1:nn
                sh(1,i)=B_(1,i*2-1);
                sh(2,i)=B_(2,i*2);
                sh(3,i)=B_(4,i*2-1);
            end
        else
            sh=zeros(2,nn);
            for i=1:nn
                sh(1,i)=B_(1,i*2-1);
                sh(2,i)=B_(2,i*2);
            end
        end

        % Stress
        for i=1:4
            sig(i,1)=Mat_state.Sigma((e-1)*4+i,1);
        end
        
        [T]=LIB.e2E(sig);
        
        % ----------------------------
        % Internal forces
        % ----------------------------
        volume=GEOMETRY.Area(e)*MAT_POINT(e).J;
        if SOLVER.AXI
            mat=[1 0 1; 0 1 0];
            Tt=mat*T;
            vol=2*pi*MAT_POINT(e).xg(1)*volume;
        else
            Tt=T(1:2,1:2);
            vol=volume;
        end
        int_forces_1=Tt*sh*vol;

        if SOLVER.UW
            sh2=sh(1:2,:);
            if SOLVER.AXI
                sh2(1,:)=sh2(1,:)+sh(3,:);
            end
            int_forces_2=sh2*Mat_state.pw(e,1)*vol;
            if SOLVER.IMPLICIT==0 && SOLVER.step0==0
                for i=1:nn
                   nod=nd(i);
                   for j=1:sp
                        Mat_state.fint(nod*df+1-sp-j,1)=...
                            Mat_state.fint(nod*df+1-sp-j,1)-int_forces_1(3-j,i);
                        Mat_state.fint(nod*df+1-j,1)=...
                            Mat_state.fint(nod*df+1-j,1)-int_forces_2(3-j,i);
                   end
                end
            else
                for i=1:nn
                   nod=nd(i);
                   for j=1:sp
                        Mat_state.fint(nod*df-1-j,1)=...
                            Mat_state.fint(nod*df-1-j,1)-int_forces_1(3-j,i)+...
                            int_forces_2(3-j,i);
                        Mat_state.fint(nod*df+1-j,1)=...
                            Mat_state.fint(nod*df+1-j,1)+int_forces_2(3-j,i);
                   end
                end
            end            
        else
            for i=1:nn
               nod=nd(i);
               for j=1:sp
                    Mat_state.fint(nod*sp+1-j,1)=Mat_state.fint(nod*sp+1-j,1)+...
                        int_forces_1(3-j,i);
               end
            end
        end
        clear sh
    end
end

    