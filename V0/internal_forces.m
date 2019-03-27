function [Mat_state]=internal_forces(Shape_function,Mat_state)

    global GEOMETRY SOLVER
    
    sig=zeros(4,1);
    Mat_state.fint(:,1) = zeros(GEOMETRY.nodes*GEOMETRY.df,1);
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;

    for e=1:GEOMETRY.elements
        
        nd = Shape_function.near{e};
        nn=length(nd);

        % Derivatives
        B_=Shape_function.B{e};
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
        
        [T]=e2E(sig);
        
        % ----------------------------
        % Internal forces
        % ----------------------------
        volume=GEOMETRY.Area(e)*Mat_state.J(e);
        if SOLVER.AXI
            mat=[1 0 1; 0 1 0];
            Tt=mat*T;
            vol=2*pi*Mat_state.xg(e,1)*volume;
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
            if SOLVER.IMPLICIT==0
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

function [E]=e2E(e)
   
    E=zeros(3,3);

    %Build matrix
    E(1,1)=e(1);
    E(2,2)=e(2);
    E(3,3)=e(3);
    E(1,2)=e(4);
    E(2,1)=e(4);
end

    