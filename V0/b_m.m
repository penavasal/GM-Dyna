
function [B]=b_m(elements,near,dp,sp,wrap,B,xg,p)

    global SOLVER
    
    if SOLVER.AXI
        for i=1:elements
            if wrap(i)==1
                r=xg(i,1);
                nb=near{i};
                sh=dp{i};
                pp=p{i};
                n=length(nb);
                b2=zeros(4,sp*n);
                for j=1:n
                    b2(1,sp*j-1)=sh(j,1);
                    b2(2,sp*j)  =sh(j,2);
                    b2(4,sp*j-1)=pp(j)/r;
                    b2(3,sp*j-1)=sh(j,2);
                    b2(3,sp*j)  =sh(j,1);
                end
                B(i)={b2};
                clear b2;
            end
        end
    else
        for i=1:elements
            if wrap(i)==1
                nb=near{i};
                sh=dp{i};
                n=length(nb);
                b=zeros(3,sp*n);
                for j=1:n
                    b(1,sp*j-1)=sh(j,1);
                    b(2,sp*j)  =sh(j,2);
                    b(3,sp*j-1)=sh(j,2);
                    b(3,sp*j)  =sh(j,1);
                end
                B(i)={b};
                clear b;
            end
        end
    end

end