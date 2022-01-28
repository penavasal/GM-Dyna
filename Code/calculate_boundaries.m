
function [boundary,i_disp,velo,matrix]=calculate_boundaries(STEP,matrix)

    global BOUNDARY
    
    ste = STEP.ste;
    BLCK= STEP.BLCK;
    %dt=STEP.dt;

    constrains  = BOUNDARY{BLCK}.constrains;
    dad         = BOUNDARY{BLCK}.dad;
    vad         = BOUNDARY{BLCK}.vad;
    b_mult      = BOUNDARY{BLCK}.b_mult; 
    TYPE        = BOUNDARY{BLCK}.Type;
       
    [dofs,bds]=size(constrains);
    boundary=zeros(dofs,1);
    i_disp=zeros(dofs,1);
    velo=zeros(dofs,1);
    
    for m=1:bds
        if TYPE(m)~=7
            t   = STEP.t;
            if t<eval(b_mult(2,m)) || t>eval(b_mult(3,m))
                continue;
            else
                if strcmp(b_mult(4,m),'VALUE')
                    val=eval(b_mult(1,m));
                elseif strcmp(b_mult(4,m),'FUNCTION')
                    val=eval(b_mult(1,m));
                else
                    disp('Error, wrong type of boundary!')
                    stop
                end 

                if ste>1
                    t=t-STEP.dt;
                    if BOUNDARY{BLCK}.Type(m)~=5
                        if t<eval(b_mult(2,m)) || t>eval(b_mult(3,m))
                            val2=0;
                        else
                            if strcmp(b_mult(4,m),'VALUE')
                                val2=eval(b_mult(1,m));
                            elseif strcmp(b_mult(4,m),'FUNCTION')
                                val2=eval(b_mult(1,m));
                            else
                                disp('Error, wrong type of boundary!')
                                stop
                            end 
                        end
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
                            %i_disp(i)=velo(i)*dt;
                        end
                    elseif constrains(i,m)==3
                        if isscalar(matrix)
                            break;
                        else
                            matrix(i,i)=matrix(i,i)+val*dad(i,m);
                            matrix(i,BOUNDARY{BLCK}.tied(i,m))=...
                                matrix(i,BOUNDARY{BLCK}.tied(i,m))-val*dad(i,m);
                        end
                    elseif constrains(i,m)==4
                        boundary(i)=4;
                    end
                end
            end
        end
    end

end