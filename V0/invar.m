

function [P2,Q2]=invar(Ss,e)

    ss=zeros(4,1);
    for i=1:4
        ss(5-i)=Ss(e*4+1-i);
    end
    
    Sc=e2E(ss);
    
    P2=(Sc(1,1)+Sc(2,2)+Sc(3,3))/3;
    s=Sc-P2*eye(3);
    Q2= s(1,1)^2 + s(2,2)^2 + s(3,3)^2 +...
       s(2,1)^2 + s(1,2)^2 + ...
       s(3,1)^2 + s(1,3)^2 + ...
       s(2,3)^2 + s(3,2)^2;
    Q2= sqrt(3/2*Q2);
end

    
function [E]=e2E(e)
   
    E=zeros(3,3);
    
    %Build matrix
    
    E(1,1)=e(1);
    E(2,2)=e(2);
    E(1,2)=e(4);
    E(2,1)=e(4);
    E(3,3)=e(3);
    
end