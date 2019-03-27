
function [xg]=update_mp(d,Shape_function,xg)
    
    global GEOMETRY
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    
    %% Update **********************
    for i=1:GEOMETRY.elements
        nd = Shape_function.near{i};
        m  = length(nd);
        sh = Shape_function.p{i};
        for t1=1:m
            for k=1:sp
                xg(i,k)=xg(i,k)...
                +sh(t1)*(d((nd(t1)-1)*df+k,1)-d((nd(t1)-1)*df+k,2));
            end
        end
    end  
end