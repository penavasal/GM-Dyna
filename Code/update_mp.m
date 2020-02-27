
function [MAT_POINT]=update_mp(d,MAT_POINT)
    
    global GEOMETRY SOLVER
    
    sp=GEOMETRY.sp;
    df=GEOMETRY.df;
    
    [xg]=LIB.S2list(MAT_POINT{1},'xg');
    %% Update **********************
    for i=1:GEOMETRY.mat_points
        nd = MAT_POINT{1}(i).near;
        m  = length(nd);
        sh = MAT_POINT{1}(i).N;
        for t1=1:m
            for k=1:sp
                xg(i,k)=xg(i,k)...
                +sh(t1)*(d((nd(t1)-1)*df+k,1)-d((nd(t1)-1)*df+k,2));
            end
        end
    end
    [MAT_POINT{1}]=LIB.list2S(MAT_POINT{1},'xg',xg);

    %% Other phases
    [phases,~]=size(SOLVER.PHASES);
    for ph=2:phases
        MAT_POINT{ph}(mp).xg=MAT_POINT{1}(mp).xg;
    end
end