
function [boundary,i_disp,velo]=calculate_boundaries(ste)

    global BOUNDARY

    constrains  = BOUNDARY.constrains;
    dad         = BOUNDARY.dad;
    vad         = BOUNDARY.vad;
    b_mult      = BOUNDARY.b_mult; 
       
    [dofs,bds]=size(constrains);
    boundary=zeros(dofs,1);
    i_disp=zeros(dofs,1);
    velo=zeros(dofs,1);
    
    for m=1:bds
        if b_mult(ste,m)=='NULL'
            continue;
        else
            val=str2double(b_mult(ste,m));
            if ste>1
                if b_mult(ste-1,m)=='NULL'
                    val2=0;
                else
                    val2=str2double(b_mult(ste-1,m));
                end
                if val2
                    val=val-val2;
                end
            end
            for i=1:dofs
                if constrains(i,m)==1
                    if boundary(i)==1
                        aux_disp=val*dad(i,m);
                        if aux_disp~=i_disp(i)
                            disp('Duplicated boundary conditions');
                            stop;
                        end
                    else
                        boundary(i)=constrains(i,m);
                        i_disp(i)=val*dad(i,m);
                    end
                elseif constrains(i,m)==2
                    if boundary(i)==1
                        aux_velo=val*vad(i,m);
                        if aux_velo~=velo(i)
                            disp('Duplicated boundary conditions');
                            stop;
                        end
                    else
                        velo(i)=val*vad(i,m);
                    end
                end
            end
        end
    end

end