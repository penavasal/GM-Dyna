classdef W_retention
    methods(Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Runge-Kutta 4 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [hy]=rk4(tinic,tfin,y0,G)

            %func='SWCC_hysteresis';
            N=500;
            
            %paso del metodo%
            hx=(tfin-tinic)/N;
            hy=0;
            z=y0;
            s=tinic;

            for i=0:N-1

               F1=W_retention.SWCC_hysteresis(s,z,hx,hy,G);
               F2=W_retention.SWCC_hysteresis(s+hx/2,z+(hx/2)*F1,hx,hy,G);
               F3=W_retention.SWCC_hysteresis(s+hx/2,z+(hx/2)*F2,hx,hy,G);
               F4=W_retention.SWCC_hysteresis(s+hx,z+hx*F3,hx,hy,G);

               z=z+hx*1/6*(F1+2*F2+2*F3+F4);   
               if z<=0.06
                   z;
               end
               s=s+hx;

            end
            hy=z-y0;
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Hysteresis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function z=SWCC_hysteresis(x,y,hx,hy,G)

            lambdad = G(1);
            yR      = G(2);
            Betad   = G(3);
            Betaw   = G(4);
            Beta1   = G(5);
            c1d     = G(6);
            c2d     = G(7);
            c3d     = G(8);
            c1w     = G(9);
            c2w     = G(10);
            c3w     = G(11);

            if (hx>0)||(hy<0)

                if hy <0
                    x1=x;
                    x=y;
                    y=x1;
                end

                yd=-lambdad*x+(1/Betad)*log(c3d+c2d*exp(c1d*x));
                Dd=y-yR;
                lambdadbar=lambdad*(1-exp(-Betad*Dd));
                Ddbar=yd-y;
                Betawbar=Betaw*sqrt(y);
                lambdabar=lambdadbar*exp(-Betawbar*Ddbar);

                if hx>0
                    z=-lambdabar;
                elseif hy<0
                    z=-1/lambdabar;
                end


            elseif (hx<0)||(hy>0)

                if hy >0
                    x1=x;
                    x=y;
                    y=x1;
                end

                yw=-lambdad*x-(1/Betaw)*log(c3w+c2w*exp(c1w*x));
                Dw=1-y;
                lambdawbar=lambdad*(1-exp(-Betaw*Dw));
                Dwbar=y-yw;
                lambdabar=lambdawbar*exp(-Beta1*Dwbar);

                if hx<0
                    z=-lambdabar;
                elseif hy>0
                    z=-1/lambdabar;
                end
                
            else
                z=0;
            end
            
            if ~isreal(z)
                z;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perm driver
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function krw=K(e,BLCK,sw)
            
            global MATERIAL GEOMETRY
    
            Material=GEOMETRY.material;
            MAT=MATERIAL(BLCK).MAT;
            
            RC=MAT{61,Material(e)};
            
            if strcmp(RC,'PEDROSO')
                lambda0=MAT{50,Material(e)};
                lambda1=MAT{51,Material(e)};
                beta=MAT{52,Material(e)};
                alfa=MAT{53,Material(e)};

                c1=beta*(lambda0-lambda1);
                c2=exp(-alfa*beta);
                c3=exp(beta*(lambda0-1))-c2*exp(c1);

                sw2=sw(e,1);

                krw=lambda0*sw2-1/beta*log(c3+c2*exp(c1*sw2));
            elseif strcmp(RC,'VG')
                
                xi=MAT{50,Material(e)};
                m=MAT{58,Material(e)};
                swr=MAT{56,Material(e)};
                
                sw2=sw(e,1);
                Swe=(sw2 - swr)/(1-swr);
                
                krw = Swe^xi * (1-( 1 - Swe^(1/m))^m )^2;
            else
                error('oh, oh');
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Q driver
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function Qp=Q(e,Mat_state,K_s,K_w,n)
            
            pw=Mat_state.pw(e,:);
            sw=Mat_state.sw(e,:);
            
            dpw=pw(1)-pw(2);
            if abs(dpw)<1e-9
                cs=0;
            else
                cs=n*(sw(1)-sw(2))/dpw;
            end
            
            if e==60
                e;
            end
            
            if cs<0
                e;
            end
            
            Qp=1/(cs + sw(1)*n/K_w + sw(1)*sw(1)*(1-n)/K_s);
%             if ~isreal(Qp)
%                 Qp;
%             end

            Mat_state.cs(e,1)=cs;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sw driver
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function SW=SW(e,Mat_state,MAT,Mat)
            
            global SOLVER GEOMETRY
    
            Material=GEOMETRY.material;
            
            RC=MAT{61,Material(e)};
            
            pw=Mat_state.pw(e,:);
            sw=Mat_state.sw(e,:);
            
            % Change to kPa
            if strcmp(SOLVER.UNITS,'MM')
                pw=-pw*1000;
            elseif strcmp(SOLVER.UNITS,'M')
                pw=-pw/1000;
            else
                error('not implemented yet')
            end
            
            if strcmp(RC,'PEDROSO')

                G=W_retention.data(e,MAT,Mat);


                if pw(1)<1 || pw(2)<1
                    SW=1;             
                else
                    tinic=log(1+pw(2));
                    tfin =log(1+pw(1));
                    [hy]=W_retention.rk4(tinic,tfin,sw(2),G);
                    SW = sw(2)+hy;

                    if ~isreal(SW) || SW>1 || SW<0
                        hy;
                    end

                    if e==60
                        e;
                    end

                end
            elseif strcmp(RC,'VG')
                if pw(1)<0
                    SW=1;
                else
                    p0=MAT{55,Material(e)};
                    alfa=MAT{54,Material(e)};
                    n=MAT{60,Material(e)};
                    m=MAT{58,Material(e)};
                    swr=MAT{56,Material(e)};

                    if alfa==0
                        if p0==0
                            error('oh, oh');
                        else
                            alfa=1/p0;
                        end
                    end

                    Swe = (1 + (alfa*pw(1))^n)^(-m);
                    SW = swr + Swe*(1-swr);
                end
                
            else 
                error('oh, oh');
            end
            

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % G data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
        function G=data(e,MAT,Mat)
            
            y0=1;
            
            lambdad=MAT{54,Mat(e)};
            xRd=MAT{55,Mat(e)};
            yR=MAT{56,Mat(e)};
            xRw=MAT{57,Mat(e)};
            Betad=MAT{58,Mat(e)};
            Betaw=MAT{59,Mat(e)};
            Beta1=MAT{60,Mat(e)};

            c1d=Betad*lambdad; c2d=exp(Betad*yR); c3d=exp(Betad*(y0+lambdad*xRd))-c2d*exp(c1d*xRd);
            c1w=-Betaw*lambdad; c2w=exp(-Betaw*y0); c3w=exp(-Betaw*lambdad*xRw)-c2w*exp(c1w*xRw);

            G(1)=lambdad;
            G(2)=yR;
            G(3)=Betad;
            G(4)=Betaw;
            G(5)=Beta1;
            G(6)=c1d;
            G(7)=c2d;
            G(8)=c3d;
            G(9)=c1w;
            G(10)=c2w;
            G(11)=c3w;
        end
    end
end

