

function element_FEM()

    
    
    [dxidx, N, dNdx, w,Dxdxi] = elem_driver(coord,xilist);


end



function [xilist]=search_xilist(coord,int_point_pos,xilist,i)
    
    [dxidx_0,~,~,~,~] = elem_driver(coord,xilist(i-1,:));
     dx=int_point_pos(i,:)-int_point_pos(i-1,:);
     xilist(i,:)=xilist(i-1,:)+dx*dxidx_0;
     
     
end   

function [xi]=search_xilist_NR(coord,int_point_pos,xi_0,coord_xi)
    
    [NNE,~]=size(coord_xi);
    
    if NNE==3
        elem=[1 2 3];
    else
        elem=[1 2 3 4];
    end

    tol = 1e-3;
    xi=xi_0;
    norm_f(1)=1000;
    norm_r_f(1)=1000;
    conver=0;
    a=1;
    i=1;
    while conver==0
        i=i+1;
         [~, N, dNdx,~,Dxdxi] = ...
             elem_driver(coord,xi);

         f=int_point_pos-N'*coord;
         df=-Dxdxi;
         %df=-dNdx'*coord;
         xi=xi_0-a*f/df;
         
         [I]=IoO(xi,coord_xi,elem);
         while I==0
             [I]=IoO(xi,coord_xi,elem);
             if I==0
                a=a/2;
                xi=xi_0-a*f/df;
             else
                 break;
             end
         end
         a=1;
        norm_f(i)=norm(f);
        norm_r_f(i)=abs(norm_f(i)-norm_f(i-1));
        
        if abs(norm_f(i))<tol || abs(norm_r_f(i))<tol/100
            conver=1;
        end
        xi_0=xi;
    end
    
    xi;

end
%
%========================== ELEMENT DRIVER ================================
%
function [dxidx, N, dNdx, w,Dxdxi] = elem_driver(coord,xilist)
%
%
%    Arguments;
%
%      xi(i,inpt)         local coord of integration point no. intpt
%      coords(i,a)        ith coord of ath node
%
%   Local variables
%      npoints            No. integration points
%      w(intpt)           weight for integration point no. intpt
%      N(a)               Shape function associated with ath node on element
%      dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
%      dNdx(a,i)          Derivative of ath shape function wrt ith global coord
%      dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
%      dxidx(i,j)         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian
%


%      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
%      nelnodes           No. nodes on the element
   [nelnodes,ncoord]=size(coord);
   [npoints,~]=size(xilist);
   dNdx = zeros(nelnodes,ncoord);
   dxdxi = zeros(ncoord,ncoord);
%
%  Set up integration points && weights    
   w = integrationweights(ncoord,nelnodes,npoints);
%
%  Loop over the integration points
%
   for intpt = 1:npoints

%     Compute shape functions && derivatives wrt local coords
%
      xi=zeros(ncoord,1);
      for i = 1:ncoord
        xi(i) = xilist(intpt,i);
      end      
      
      N = shapefunctions(nelnodes,ncoord,xi);
      dNdxi = shapefunctionderivs(nelnodes,ncoord,xi);

      
%
%     Compute the jacobian matrix && its determinant
%
      for i = 1:ncoord
        for j = 1:ncoord
          dxdxi(i,j) = 0.;
          for a = 1:nelnodes
            dxdxi(i,j) = dxdxi(i,j) + coord(a,i)*dNdxi(a,j);
          end
        end
      end
      
      Dxdxi=transpose(dxdxi);
      dxidx = inv(transpose(dxdxi));
      dt = det(dxdxi);
      
      %
%     Convert shape function derivatives:derivatives wrt global coords
%
      for a = 1:nelnodes
        for i = 1:ncoord
          dNdx(a,i) = 0.;
          for j = 1:ncoord
            dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
          end
        end
      end

       w(intpt)=w(intpt)*dt;

   end

end

%
%================= INTEGRATION WEIGHTS ==================================
%
%   Defines integration weights w_i
%
function w = integrationweights(ncoord,nelnodes,npoints)

   w = zeros(npoints,1);

%
%  1D elements
%
   if (ncoord == 1) 
     if (npoints == 1)
       w(1) = 2.;
     elseif (npoints == 2) 
       w = [1.,1.];
     elseif (npoints == 3) 
       w = [0.555555555,0.888888888,0.555555555];
     end
%
%  2D elements
%
   elseif (ncoord == 2) 
%
%    Triangular element
%
     if ( nelnodes == 3 || nelnodes == 6 ) 
       if (npoints == 1) 
         w(1) = 0.5;
       elseif (npoints == 3) 
         w(1) = 1./6.;
         w(2) = 1./6.;
         w(3) = 1./6.;
       elseif (npoints == 4) 
         w = [-27./96.,25./96.,25/96.,25/96.];
       end
%
%    Rectangular element
%                  
     elseif ( nelnodes==4 || nelnodes==8 ) 

       if (npoints == 1) 
         w(1) = 4.;
       elseif (npoints == 4) 
         w = [1.,1.,1.,1.];
       elseif (npoints == 9 ) 
         w1D = [0.555555555,0.888888888,0.55555555555];
         for j = 1:3
           for i = 1:3
             n = 3*(j-1)+i;
             w(n) = w1D(i)*w1D(j);
           end
         end    
       end
     end 

   elseif (ncoord == 3) 
%
%  3D elements
%
     if (nelnodes == 4 || nelnodes==10 ) 
       if (npoints == 1) 
         w(1) = 1./6.;
       elseif (npoints == 4) 
         w = [1./24.,1./24.,1./24.,1./24.];
       end
     elseif ( nelnodes==8 || nelnodes==20 ) 
       if (npoints == 1) 
         w(1) = 8.;
       elseif (npoints == 8) 
         w = [1.,1.,1.,1.,1.,1.,1.,1.];
       elseif (npoints == 27) 
         w1D = [0.555555555,0.888888888,0.55555555555];
         for k = 1:3
           for j = 1:3
             for i = 1:3
               n = 9*(k-1)+3*(j-1)+i;
               w(n) = w1D(i)*w1D(j)*w1D(k);
             end
           end    
         end
       end
     end
   end
end

%
%================= SHAPE FUNCTIONS ==================================
%
%        Calculates shape functions for various element types
%
function N = shapefunctions(nelnodes,ncoord,xi)
 

   N = zeros(nelnodes,1);
%
%  1D elements
%
  if (ncoord == 1) 
    if (nelnodes==2) 
      N(1) = 0.5*(1.+xi(1));
      N(2) = 0.5*(1.-xi(1));
    elseif (nelnodes == 3) 
      N(1) = -0.5*xi(1)*(1.-xi(1));
      N(2) =  0.5*xi(1)*(1.+xi(1));
      N(3) = (1.-xi(1))*(1.+xi(1));
    end
%
%  2D elements
%
   elseif (ncoord == 2) 
%
%    Triangular element
%
     if ( nelnodes == 3 ) 
       N(1) = xi(1);
       N(2) = xi(2);
       N(3) = 1.-xi(1)-xi(2);               
     elseif ( nelnodes == 6 ) 
       xi3 = 1.-xi(1)-xi(2);
       N(1) = (2.*xi(1)-1.)*xi(1);
       N(2) = (2.*xi(2)-1.)*xi(2);
       N(3) = (2.*xi3-1.)*xi3;
       N(4) = 4.*xi(1)*xi(2);
       N(5) = 4.*xi(2)*xi3;
       N(6) = 4.*xi3*xi(1);
%
%    Rectangular element
%                  
     elseif ( nelnodes == 4 ) 
       N(1) = 0.25*(1.-xi(1))*(1.-xi(2));
       N(2) = 0.25*(1.+xi(1))*(1.-xi(2));
       N(3) = 0.25*(1.+xi(1))*(1.+xi(2));
       N(4) = 0.25*(1.-xi(1))*(1.+xi(2));
     elseif (nelnodes == 8) 
       N(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
       N(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
       N(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
       N(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
       N(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
       N(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
       N(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
       N(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
     end
%
   elseif (ncoord==3) 

     if (nelnodes == 4) 
       N(1) = xi(1);
       N(2) = xi(2);
       N(3) = xi(3);
       N(4) = 1.-xi(1)-xi(2)-xi(3);
     elseif (nelnodes == 10) 
       xi4 = 1.-xi(1)-xi(2)-xi(3);
       N(1) = (2.*xi(1)-1.)*xi(1);
       N(2) = (2.*xi(2)-1.)*xi(2);
       N(3) = (2.*xi(3)-1.)*xi(3);
       N(4) = (2.*xi4-1.)*xi4;
       N(5) = 4.*xi(1)*xi(2);
       N(6) = 4.*xi(2)*xi(3);
       N(7) = 4.*xi(3)*xi(1);
       N(8) = 4.*xi(1)*xi4;
       N(9) = 4.*xi(2)*xi4;
       N(10) = 4.*xi(3)*xi4;
     elseif (nelnodes == 8) 
       N(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))/8.;
       N(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))/8.;
       N(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))/8.;
       N(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))/8.;
       N(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))/8.;
       N(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))/8.;
       N(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))/8.;
       N(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))/8.;
     elseif (nelnodes == 20) 
       N(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)/8.;
       N(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)/8.;
       N(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)/8.;
       N(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)/8.;
       N(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)/8.;
       N(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)/8.;
       N(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)/8.;
       N(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)/8.;
       N(9)  = (1.-xi(1)^2)*(1.-xi(2))*(1.-xi(3))/4.;
       N(10) = (1.+xi(1))*(1.-xi(2)^2)*(1.-xi(3))/4.;
       N(11) = (1.-xi(1)^2)*(1.+xi(2))*(1.-xi(3))/4.;
       N(12) = (1.-xi(1))*(1.-xi(2)^2)*(1.-xi(3))/4.;
       N(13) = (1.-xi(1)^2)*(1.-xi(2))*(1.+xi(3))/4.;
       N(14) = (1.+xi(1))*(1.-xi(2)^2)*(1.+xi(3))/4.;
       N(15) = (1.-xi(1)^2)*(1.+xi(2))*(1.+xi(3))/4.;
       N(16) = (1.-xi(1))*(1.-xi(2)^2)*(1.+xi(3))/4.;
       N(17) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3)^2)/4.;
       N(18) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3)^2)/4.;
       N(19) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3)^2)/4.;
       N(20) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3)^2)/4.;
     end
   end

end

%
%====================== INTEGRATION POINTS ==================================
%
%   Defines positions of integration points
%
function [xi,coord_xi] = integrationpoints(ncoord,nelnodes,npoints)

   xi = zeros(ncoord,npoints);
   coord_xi = zeros(nelnodes,ncoord);
   
%
%  1D elements
%
   if (ncoord == 1) 
     if (npoints==1) 
       xi(1,1) = 0.;
     elseif (npoints == 2) 
       xi(1,1) = -0.5773502692;
       xi(1,2) = -xi(1,1);
     elseif (npoints == 3) 
       xi(1,1) = -0.7745966692;
       xi(1,2) = 0.0;
       xi(1,3) = -xi(1,1);
     end
%
%  2D elements
%
   elseif (ncoord == 2) 
%
%    Triangular element
%
     if ( nelnodes == 3 || nelnodes == 6 ) 
       if (npoints == 1) 
         xi(1,1) = 1./3.;
         xi(2,1) = 1./3.;
         
         coord_xi=[0 1;
                    1 0;
                    0 0];
                
       elseif (npoints == 3) 
         xi(1,1) = 0.6;
         xi(2,1) = 0.2;
         xi(1,2) = 0.2;
         xi(2,2) = 0.6;
         xi(1,3) = 0.2;
         xi(2,3) = 0.2;
       elseif (npoints == 4) 
         xi(1,1) = 1./3.;
         xi(2,1) = 1./3.;
         xi(1,2) = 0.6;
         xi(2,2) = 0.2;
         xi(1,3) = 0.2;
         xi(2,3) = 0.6;
         xi(1,4) = 0.2;
         xi(2,4) = 0.2;
       elseif (npoints == 7)
         xi(1,1) = 1./3.0;
  	     xi(2,1) = xi(1,1);
	     xi(1,2) = 0.0597158717;
	     xi(2,2) = 0.4701420641;
         xi(1,3) = xi(2,2);
	     xi(2,3) = xi(1,2);
	     xi(1,4) = xi(2,2);
	     xi(2,4) = xi(2,2);
	     xi(1,5) = 0.7974269853;
	     xi(2,5) = 0.1012865073;
	   	 xi(1,6) = xi(2,5);
	     xi(2,6) = xi(1,5);
	     xi(1,7) = xi(2,5);
	     xi(2,7) = xi(2,5);
       end
%
%    Rectangular element
%                  
     elseif ( nelnodes==4 || nelnodes==8 ) 

       if (npoints == 1) 
         xi(1,1) = 0.;
         xi(2,1) = 0.;
         
         coord_xi=[-1 -1;
                    1 -1;
                    1 1;
                    -1 1];
                
       elseif (npoints == 4) 
         xi(1,1) = -0.5773502692;
         xi(2,1) = xi(1,1);
         xi(1,2) = -xi(1,1);
         xi(2,2) = xi(1,1);
         xi(1,3) = xi(1,1);
         xi(2,3) = -xi(1,1);
         xi(1,4) = -xi(1,1);
         xi(2,4) = -xi(1,1);
       elseif (npoints == 9) 
         xi(1,1) = -0.7745966692;
         xi(2,1) = xi(1,1);
         xi(1,2) = 0.0;
         xi(2,2) = xi(1,1);
         xi(1,3) = -xi(1,1);
         xi(2,3) = xi(1,1);
         xi(1,4) = xi(1,1);
         xi(2,4) = 0.0;
         xi(1,5) = 0.0;
         xi(2,5) = 0.0;
         xi(1,6) = -xi(1,1);
         xi(2,6) = 0.0;
         xi(1,7) = xi(1,1);
         xi(2,7) = -xi(1,1);
         xi(1,8) = 0.;
         xi(2,8) = -xi(1,1);
         xi(1,9) = -xi(1,1);
         xi(2,9) = -xi(1,1);
       end
     end
%
%   3D elements
%
   elseif (ncoord == 3) 
%
%  3D elements
%
     if (nelnodes == 4 || nelnodes==10 ) 
       if (npoints == 1) 
         xi(1,1) = 0.25;
         xi(2,1) = 0.25;
         xi(3,1) = 0.25;
       elseif (npoints == 4) 
         xi(1,1) = 0.58541020;
         xi(2,1) = 0.13819660;
         xi(3,1) = xi(2,1);
         xi(1,2) = xi(2,1);
         xi(2,2) = xi(1,1);
         xi(3,2) = xi(2,1);
         xi(1,3) = xi(2,1);
         xi(2,3) = xi(2,1);
         xi(3,3) = xi(1,1);
         xi(1,4) = xi(2,1);
         xi(2,4) = xi(2,1);
         xi(3,4) = xi(2,1);
        elseif (npoints == 5)
         xi(1,1) = 0.25;
         xi(2,1) = 0.25;
         xi(3,1) = 0.25;
         xi(1,2) = 0.5;
         xi(2,2) = 1./6.;
         xi(3,2) = 1./6.;
         xi(1,3) = 1./6.;
         xi(2,3) = 0.5;
         xi(3,3) = 1./6.;
         xi(1,4) = 1./6.;
         xi(2,4) = 1./6.;
         xi(3,4) = 0.5;
         xi(1,5) = 1./6.;
         xi(2,5) = 1./6.;
         xi(3,5) = 1./6.;
       end
     elseif ( nelnodes==8 || nelnodes==20 ) 
       if (npoints == 1) 
         xi(1,1) = 0.;
         xi(2,1) = 0.;
         xi(3,1) = 0.;
       elseif (npoints == 8) 
         x1D = [-0.5773502692,0.5773502692];
         for k = 1:2
           for j = 1:2 
             for i = 1:2
               n = 4*(k-1) + 2*(j-1) + i;
               xi(1,n) = x1D(i);
               xi(2,n) = x1D(j);
               xi(3,n) = x1D(k);
             end
           end
         end
       elseif (npoints == 27) 
         x1D = [-0.7745966692,0.,0.7745966692];
         for k = 1:3
           for j = 1:3
             for i = 1:3
               n = 9*(k-1) + 3*(j-1) + i;
               xi(1,n) = x1D(i);
               xi(2,n) = x1D(j);
               xi(3,n) = x1D(k);
             end
           end
         end
       end
     end
   end
   
   xi=transpose(xi);
end

%
%================= SHAPE FUNCTION DERIVATIVES ======================
%
function dNdxi = shapefunctionderivs(nelnodes,ncoord,xi)

  dNdxi = zeros(nelnodes,ncoord);
%
% 1D elements
%
  if (ncoord == 1) 
    if (nelnodes==2) 
      dNdxi(1,1) = 0.5;
      dNdxi(2,1) = -0.5;
    elseif (nelnodes == 3) 
      dNdxi(1,1) = -0.5+xi(1);
      dNdxi(2,1) =  0.5+xi(1);
      dNdxi(3,1) = -2.*xi(1);
    end
%
%  2D elements
%
   elseif (ncoord == 2) 
%
%    Triangular element
%
     if ( nelnodes == 3 ) 
       dNdxi(1,1) = 1.;
       dNdxi(2,2) = 1.;
       dNdxi(3,1) = -1.;
       dNdxi(3,2) = -1.;               
     elseif ( nelnodes == 6 ) 
       xi3 = 1.-xi(1)-xi(2);
       dNdxi(1,1) = 4.*xi(1)-1.;
       dNdxi(2,2) = 4.*xi(2)-1.;
       dNdxi(3,1) = -(4.*xi3-1.);
       dNdxi(3,2) = -(4.*xi3-1.);
       dNdxi(4,1) = 4.*xi(2);
       dNdxi(4,2) = 4.*xi(1);
       dNdxi(5,1) = -4.*xi(2);
       dNdxi(5,2) = 4.*xi3 - 4.*xi(2);
       dNdxi(6,1) = 4.*xi3 - 4.*xi(1);
       dNdxi(6,2) = -4.*xi(1);
%
%    Rectangular element
%                  
     elseif ( nelnodes == 4 ) 
       dNdxi(1,1) = -0.25*(1.-xi(2));
       dNdxi(1,2) = -0.25*(1.-xi(1));
       dNdxi(2,1) = 0.25*(1.-xi(2));
       dNdxi(2,2) = -0.25*(1.+xi(1));
       dNdxi(3,1) = 0.25*(1.+xi(2));
       dNdxi(3,2) = 0.25*(1.+xi(1));
       dNdxi(4,1) = -0.25*(1.+xi(2));
       dNdxi(4,2) = 0.25*(1.-xi(1));
     elseif (nelnodes == 8) 
       dNdxi(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
       dNdxi(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
       dNdxi(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
       dNdxi(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
       dNdxi(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
       dNdxi(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
       dNdxi(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
       dNdxi(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
       dNdxi(5,1) = -xi(1)*(1.-xi(2));
       dNdxi(5,2) = -0.5*(1.-xi(1)*xi(1));
       dNdxi(6,1) = 0.5*(1.-xi(2)*xi(2));
       dNdxi(6,2) = -(1.+xi(1))*xi(2);
       dNdxi(7,1) = -xi(1)*(1.+xi(2));
       dNdxi(7,2) = 0.5*(1.-xi(1)*xi(1));
       dNdxi(8,1) = -0.5*(1.-xi(2)*xi(2));
       dNdxi(8,2) = -(1.-xi(1))*xi(2);
      end
%
%    3D elements
%
   elseif (ncoord==3) 

     if (nelnodes == 4) 
       dNdxi(1,1) = 1.;
       dNdxi(2,2) = 1.;
       dNdxi(3,3) = 1.;
       dNdxi(4,1) = -1.;
       dNdxi(4,2) = -1.;
       dNdxi(4,3) = -1.;
     elseif (nelnodes == 10) 
       xi4 = 1.-xi(1)-xi(2)-xi(3);
       dNdxi(1,1) = (4.*xi(1)-1.);
       dNdxi(2,2) = (4.*xi(2)-1.);
       dNdxi(3,3) = (4.*xi(3)-1.);
       dNdxi(4,1) = -(4.*xi4-1.);
       dNdxi(4,2) = -(4.*xi4-1.);
       dNdxi(4,3) = -(4.*xi4-1.);
       dNdxi(5,1) = 4.*xi(2);
       dNdxi(5,2) = 4.*xi(1);
       dNdxi(6,2) = 4.*xi(3);
       dNdxi(6,3) = 4.*xi(2);
       dNdxi(7,1) = 4.*xi(3);
       dNdxi(7,3) = 4.*xi(1); 
       dNdxi(8,1) = 4.*(xi4-xi(1));
       dNdxi(8,2) = -4.*xi(1);
       dNdxi(8,3) = -4.*xi(1);
       dNdxi(9,1) = -4.*xi(2);
       dNdxi(9,2) = 4.*(xi4-xi(2));
       dNdxi(9,3) = -4.*xi(2);
       dNdxi(10,1) = -4.*xi(3)*xi4;
       dNdxi(10,2) = -4.*xi(3);
       dNdxi(10,3) = 4.*(xi4-xi(3));
     elseif (nelnodes == 8) 
       dNdxi(1,1) = -(1.-xi(2))*(1.-xi(3))/8.;
       dNdxi(1,2) = -(1.-xi(1))*(1.-xi(3))/8.;
       dNdxi(1,3) = -(1.-xi(1))*(1.-xi(2))/8.;
       dNdxi(2,1) = (1.-xi(2))*(1.-xi(3))/8.;
       dNdxi(2,2) = -(1.+xi(1))*(1.-xi(3))/8.;
       dNdxi(2,3) = -(1.+xi(1))*(1.-xi(2))/8.;
       dNdxi(3,1) = (1.+xi(2))*(1.-xi(3))/8.;
       dNdxi(3,2) = (1.+xi(1))*(1.-xi(3))/8.;
       dNdxi(3,3) = -(1.+xi(1))*(1.+xi(2))/8.;
       dNdxi(4,1) = -(1.+xi(2))*(1.-xi(3))/8.;
       dNdxi(4,2) = (1.-xi(1))*(1.-xi(3))/8.;
       dNdxi(4,3) = -(1.-xi(1))*(1.+xi(2))/8.;
       dNdxi(5,1) = -(1.-xi(2))*(1.+xi(3))/8.;
       dNdxi(5,2) = -(1.-xi(1))*(1.+xi(3))/8.;
       dNdxi(5,3) = (1.-xi(1))*(1.-xi(2))/8.;
       dNdxi(6,1) = (1.-xi(2))*(1.+xi(3))/8.;
       dNdxi(6,2) = -(1.+xi(1))*(1.+xi(3))/8.;
       dNdxi(6,3) = (1.+xi(1))*(1.-xi(2))/8.;
       dNdxi(7,1) = (1.+xi(2))*(1.+xi(3))/8.;
       dNdxi(7,2) = (1.+xi(1))*(1.+xi(3))/8.;
       dNdxi(7,3) = (1.+xi(1))*(1.+xi(2))/8.;
       dNdxi(8,1) = -(1.+xi(2))*(1.+xi(3))/8.;
       dNdxi(8,2) = (1.-xi(1))*(1.+xi(3))/8.;
       dNdxi(8,3) = (1.-xi(1))*(1.+xi(2))/8.;
     elseif (nelnodes == 20) 
       dNdxi(1,1) = (-(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;
       dNdxi(1,2) = (-(1.-xi(1))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;
       dNdxi(1,3) = (-(1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)-xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;

       dNdxi(2,1) = ((1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)+(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;
       dNdxi(2,2) = (-(1.+xi(1))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)-(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;
       dNdxi(2,3) = (-(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-xi(3)-2.)-(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;

       dNdxi(3,1) = ((1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;
       dNdxi(3,2) = ((1.+xi(1))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;
       dNdxi(3,3) = (-(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-xi(3)-2.)-(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;

       dNdxi(4,1) = (-(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)-(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;
       dNdxi(4,2) = ((1.-xi(1))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)+(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;
       dNdxi(4,3) = (-(1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)-xi(3)-2.)-(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;
       dNdxi(5,1) = (-(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(5,2) = (-(1.-xi(1))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(5,3) = ((1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)+xi(3)-2.)+(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(6,1) = ((1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)+(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(6,2) = (-(1.+xi(1))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)-(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(6,3) = ((1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)+xi(3)-2.)+(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(7,1) = ((1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(7,2) = ((1.+xi(1))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(7,3) = ((1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)+xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(8,1) = (-(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)-(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(8,2) = ((1.-xi(1))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)+(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(8,3) = ((1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)+xi(3)-2.)+(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(9,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.;
       dNdxi(9,2)  = -(1.-xi(1)^2)*(1.-xi(3))/4.;
       dNdxi(9,3)  = -(1.-xi(1)^2)*(1.-xi(2))/4.;
       dNdxi(10,1)  = (1.-xi(2)^2)*(1.-xi(3))/4.;
       dNdxi(10,2)  = -2.*xi(2)*(1.+xi(1))*(1.-xi(3))/4.;
       dNdxi(10,3)  = -(1.-xi(2)^2)*(1.+xi(1))/4.;
       dNdxi(11,1)  = -2.*xi(1)*(1.+xi(2))*(1.-xi(3))/4.;
       dNdxi(11,2)  = (1.-xi(1)^2)*(1.-xi(3))/4.;
       dNdxi(11,3)  = -(1.-xi(1)^2)*(1.+xi(2))/4.;
       dNdxi(12,1)  = -(1.-xi(2)^2)*(1.-xi(3))/4.;
       dNdxi(12,2)  = -2.*xi(2)*(1.-xi(1))*(1.-xi(3))/4.;
       dNdxi(12,3)  = -(1.-xi(2)^2)*(1.-xi(1))/4.;
       dNdxi(13,1)  = -2.*xi(1)*(1.-xi(2))*(1.+xi(3))/4.;
       dNdxi(13,2)  = -(1.-xi(1)^2)*(1.+xi(3))/4.;
       dNdxi(13,3)  = (1.-xi(1)^2)*(1.-xi(2))/4.;
       dNdxi(14,1)  = (1.-xi(2)^2)*(1.+xi(3))/4.;
       dNdxi(14,2)  = -2.*xi(2)*(1.+xi(1))*(1.+xi(3))/4.;
       dNdxi(14,3)  = (1.-xi(2)^2)*(1.+xi(1))/4.;
       dNdxi(15,1)  = -2.*xi(1)*(1.+xi(2))*(1.+xi(3))/4.;
       dNdxi(15,2)  =  (1.-xi(1)^2)*(1.+xi(3))/4.;
       dNdxi(15,3)  = (1.-xi(1)^2)*(1.+xi(2))/4.;
       dNdxi(16,1)  = -(1.-xi(2)^2)*(1.+xi(3))/4.;
       dNdxi(16,2)  = -2.*xi(2)*(1.-xi(1))*(1.+xi(3))/4.;
       dNdxi(16,3)  = (1.-xi(2)^2)*(1.-xi(1))/4.;
       dNdxi(17,1) = -(1.-xi(2))*(1.-xi(3)^2)/4.;
       dNdxi(17,2) = -(1.-xi(1))*(1.-xi(3)^2)/4.;
       dNdxi(17,3) = -xi(3)*(1.-xi(1))*(1.-xi(2))/2.;
       dNdxi(18,1) = (1.-xi(2))*(1.-xi(3)^2)/4.;
       dNdxi(18,2) = -(1.+xi(1))*(1.-xi(3)^2)/4.;
       dNdxi(18,3) = -xi(3)*(1.+xi(1))*(1.-xi(2))/2.;
       dNdxi(19,1) = (1.+xi(2))*(1.-xi(3)^2)/4.;
       dNdxi(19,2) = (1.+xi(1))*(1.-xi(3)^2)/4.;
       dNdxi(19,3) = -xi(3)*(1.+xi(1))*(1.+xi(2))/2.;
       dNdxi(20,1) = -(1.+xi(2))*(1.-xi(3)^2)/4.;
       dNdxi(20,2) = (1.-xi(1))*(1.-xi(3)^2)/4.;
       dNdxi(20,3) = -xi(3)*(1.-xi(1))*(1.+xi(2))/2.;
     end
  end
end

