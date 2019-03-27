function [near]=make_near(x_a,x_sample,range,wrap,near,Mat,Mat_nds)
    % Naive construction of the neighbor list
    [n_a,ndim]=size(x_a);
    [n_sample,~]=size(x_sample);

    for i=1:n_sample
        if wrap(i)==1
              near_=near{i};
              mm=length(near_);
              clear near(i)
              if (ndim==1)
                x=x_sample(i);
              else
                x=x_sample(i,:);
              end
              %Naive find in range
              dist=zeros(n_a,1);
              for id=1:ndim
                dist=dist + (x_a(:,id)-x(id)).^2;
              end
              dist=sqrt(dist);
              
              for j=1:n_a
                  %if dist(j)<range(i)
                  if (Mat(i)==Mat_nds(j)) && (dist(j)<range(i))
                      t=0;
                      for k=1:mm
                          if near_(k)==j
                              t=1;
                          end
                      end
                      if t==0;
                          mm=mm+1;
                          near_(mm)=j;
                      end
                  end
              end
            near(i)={near_};
        end
    end
    clear dist
end
