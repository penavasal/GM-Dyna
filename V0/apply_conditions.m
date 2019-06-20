
function [InvK,GT]=apply_conditions(type,ste,matrix,GT)

    global GEOMETRY
    
    %%%% ONLY K MATRIX
    if type==0
        I=eye(GEOMETRY.df*GEOMETRY.nodes);
        [boundary,~,~]=calculate_boundaries(ste);
        for i=1:GEOMETRY.nodes*GEOMETRY.df
            if (boundary(i)~=0)
                for j=1:GEOMETRY.nodes*GEOMETRY.df
                    matrix(i,j)=0;
                    matrix(j,i)=0;
                end
                matrix(i,i)=1;
            end
        end
        if (abs(rcond(matrix)))<1e-14
            %fprintf('fail with K matrix in ste %i \n',ste);
        end
        InvK=matrix\I;
    
    %%%% ONLY RESIDUUM VECTOR (FIRST ITER of NEWTON-RAPHSON)
    elseif type==1
    
        [boundary,i_disp,~]=calculate_boundaries(ste);

        for i=1:GEOMETRY.nodes*GEOMETRY.df
            if (boundary(i)~=0)               
                GT(i)=i_disp(i);
            else
                for j=1:GEOMETRY.nodes*GEOMETRY.df
                    if (boundary(j)~=0)
                        GT(i)=GT(i)-matrix(i,j)*i_disp(j);
                    end
                end                     
            end
        end
        InvK=0;
    
    %%%% ONLY RESIDUUM VECTOR (SECOND and MORE ITERS of NEWTON-RAPHSON)
    elseif type==2
         [boundary,~,~]=calculate_boundaries(ste);
         for i=1:GEOMETRY.nodes*GEOMETRY.df
            if (boundary(i)~=0)               
                GT(i)=0;
            else
                for j=1:GEOMETRY.nodes*GEOMETRY.df
                    if (boundary(j)~=0)
                        GT(i)=GT(i);
                    end
                end                     
            end
         end  
         InvK=0;
         
    %%%% K and GT MATRIX for the INITIAL STEP    
    elseif type==3
        I=eye(GEOMETRY.df*GEOMETRY.nodes);
        [boundary,~,~]=calculate_boundaries(1);
        for i=1:GEOMETRY.nodes*GEOMETRY.df
            if (boundary(i)~=0)
                GT(i)=0;
                for j=1:GEOMETRY.nodes*GEOMETRY.df
                    matrix(i,j)=0;
                    matrix(j,i)=0;
                end
                matrix(i,i)=1;
            end
        end
        if (abs(rcond(matrix)))<1e-14
            fprintf('fail with K matrix in ste %i \n',ste);
        end
        InvK=matrix\I;
        
    end
    
end