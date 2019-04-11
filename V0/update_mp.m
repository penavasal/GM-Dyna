
function [MAT_POINT]=update_mp(d,MAT_POINT)
    
    global GEOMETRY
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    
    [xg]=AUX.S2list(MAT_POINT,'xg');
    %% Update **********************
    for i=1:GEOMETRY.mat_points
        nd = MAT_POINT(i).near;
        m  = length(nd);
        sh = MAT_POINT(i).N;
        for t1=1:m
            for k=1:sp
                xg(i,k)=xg(i,k)...
                +sh(t1)*(d((nd(t1)-1)*df+k,1)-d((nd(t1)-1)*df+k,2));
            end
        end
    end
    [MAT_POINT]=AUX.list2S(MAT_POINT,'xg',xg);
end