function [OUTPUT]=read_output

% File: read_output
%   Read the desired output from output.txt
%
% Date:
%   Version 1.0   29.05.2018
    global GEOMETRY SOLVER
    
    %GEOM
    L=max(GEOMETRY.x_0(:,1));
    H=max(GEOMETRY.x_0(:,2));
    
    L0=min(GEOMETRY.x_0(:,1));
    H0=min(GEOMETRY.x_0(:,2));
    
    l0=GEOMETRY.df*GEOMETRY.nodes;
    

    OUTPUT.title='OUTPUT variables';
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Read output.txt
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen('output.txt', 'rt'); % opción rt para abrir en modo texto
    formato = '%s %s %s %s'; % formato de cada línea 
    data = textscan(fid, formato, 'HeaderLines', 1);
    a = data{1};
    b = data{2};
    c = data{3};
    d = data{4};
    
    name='';
    
    [l,~] = size(a);
    t=0;
    os=0;
    while (t<l)
        t=t+1;
        s1=a{t};
        switch s1
            case '//'
                continue
            case 'OUTPUT_NAME'
                name=b{t}; 
                continue
            case 'OUTPUTS'
                os = str2double(b{t});
                type=zeros(os,2);
                RANGE=zeros(GEOMETRY.sp,2*os); % Range of loads
                VECTOR=zeros(GEOMETRY.sp,os); % Directions of loads
                continue
            case 'OUTPUT'
                M=str2double(b{t});
                if M>os
                    break;
                else
                    if strcmp(c{t},'LOAD')
                        type(M,1)=1;
                    elseif strcmp(c{t},'BOUNDARY')
                        type(M,1)=2;
                    elseif strcmp(c{t},'REACTION')
                        type(M,1)=3;
                    else
                        disp('Error, type of output not implented yet!')
                        stop
                    end
                    
                    while(t<l)
                        t=t+1;
                        s2=a{t};
                        switch s2
                            case '//'
                                continue
                            case 'OUTPUT'
                                t=t-1;
                                break
                            case 'ASSOCIATED'
                                type(M,2)=str2double(b{t});
                                continue
                            case 'X_RANGE'
                                cc = str2double(b{t});
                                if isnan(cc)
                                    if strcmp(b{t},'FULL')
                                        RANGE(1,M*2-1)=L;
                                    elseif strcmp(b{t},'INI')
                                        RANGE(1,M*2-1)=L0;
                                    else
                                        disp('Error, wrong load X range!')
                                    end
                                else
                                    RANGE(1,M*2-1)=cc;
                                end
                                cc = str2double(c{t});
                                if isnan(cc)
                                    len=strlength(c{t});
                                    if len==0 
                                        RANGE(1,M*2)=RANGE(1,M*2-1);
                                    elseif strcmp(c{t},'FULL')
                                        RANGE(1,M*2)=L;
                                    elseif strcmp(c{t},'INI')
                                        RANGE(1,M*2)=L0;
                                    else
                                        disp('Error, wrong load X range!')
                                        stop
                                    end
                                else
                                    RANGE(1,M*2)=cc;
                                end
                                continue
                            case 'Y_RANGE'
                                cc = str2double(b{t});
                                if isnan(cc)
                                    if strcmp(b{t},'FULL')
                                        RANGE(2,M*2-1)=H;
                                    elseif strcmp(b{t},'INI')
                                        RANGE(2,M*2-1)=H0;
                                    else
                                        disp('Error, wrong load Y range!')
                                    end
                                else
                                    RANGE(2,M*2-1)=cc;
                                end
                                cc = str2double(c{t});
                                if isnan(cc)
                                    len=strlength(c{t});
                                    if len==0 
                                        RANGE(2,M*2)=RANGE(2,M*2-1);
                                    elseif strcmp(c{t},'FULL')
                                        RANGE(2,M*2)=H;
                                    elseif strcmp(c{t},'INI')
                                        RANGE(2,M*2)=H0;
                                    else
                                        disp('Error, wrong load Y range!')
                                        stop
                                    end
                                else
                                    RANGE(2,M*2)=cc;
                                end
                                continue
                            case 'VECTOR'
                                VECTOR(1,M)=str2double(b{t});
                                VECTOR(2,M)=str2double(c{t});
                                if GEOMETRY.sp==3
                                    VECTOR(3,M)=str2double(d{t});
                                end
                                continue
                            otherwise
                            disp('Error, type of load not implented yet!')
                            stop
                        end
                    end
                end
                continue
            otherwise
                fprintf('Error, unrecognized parameter: %s !!\n',s1)
                stop
        end
    end
    
    fclose(fid); 
   
    nds=zeros(l0,os);
    OUTPUT.number=os;
    OUTPUT.type=zeros(os,2);
    for m=1:os
        if type(m,1)==1
            OUTPUT.type(m,1)=1;
            OUTPUT.type(m,2)=type(m,2);
            %[nds(:,m)]=search_l(type(m,2));
        elseif type(m,1)==2
            OUTPUT.type(m)=2;
            [nds(:,m)]=search_b(type(m,2));
        elseif type(m,1)==3
            OUTPUT.type(m)=2;
            [nds(:,m)]=localization(RANGE(:,2*m-1:2*m),VECTOR(:,m));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save OUTPUT_param
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OUTPUT.ref_list=nds;
    OUTPUT.list=zeros(SOLVER.dim,os);
    OUTPUT.inst=zeros(os,1);

    if strcmp(name,'')
       OUTPUT.name='FILE.mat';
    else
       OUTPUT.name=strcat(name,'.mat');
    end

    A=exist(OUTPUT.name,'file');
   
    j=0;
    while A==2
        j=j+1;
        OUTPUT.name=strcat(name,'_',int2str(j),'.mat');
        A=exist(OUTPUT.name,'file');
    end

end

function [dofs]=localization(R,dir)
    global GEOMETRY
   
    nds=zeros(GEOMETRY.nodes,1);
    for nodo=1:GEOMETRY.nodes
        tol=GEOMETRY.h_nds(nodo,1)/5;
        if (GEOMETRY.x_0(nodo,2)>=R(2,1)-tol) ...
                && (GEOMETRY.x_0(nodo,2)<=R(2,2)+tol) 
            if (GEOMETRY.x_0(nodo,1)>=R(1,1)-tol) ...
                    && (GEOMETRY.x_0(nodo,1)<=R(1,2)+tol)
                nds(nodo)=1;
            end
        end
    end
    
    dofs=zeros(GEOMETRY.nodes*GEOMETRY.df,1);
    
    for i=1:GEOMETRY.nodes
        for k=1:GEOMETRY.sp
            if nds(i)==1 && dir(k)
                dofs((i-1)*GEOMETRY.df+k)=1;
            end
        end
    end

end

function [nds]=search_b(t)
    global BOUNDARY
    [a,~]=size(BOUNDARY.constrains);
    nds=zeros(a,1);
    for i=1:a
        if BOUNDARY.constrains(i,t)
            nds(i)=1;
        end
    end
end

% function [nds]=search_l(t)
%     global ext_forces_s ext_forces_w sp UW
%     
%     [a,~]=size(ext_forces_s);
%     b=a/sp;
%     
%     if UW==0
%         d=sp;
%         nds=zeros(a,1);
%     elseif UW==1
%         d=2*sp;
%         nds=zeros(2*a,1);
%     elseif UW==2
%         d=sp*2+1;
%         nds=zeros(d*a/sp,1);
%     end
%     
%     for i=1:b
%         for k=1:sp
%             if ext_forces_s(i*sp+1-k,t)
%                 nds(i*d+1-sp-k)=1;
%             elseif UW==1 && ext_forces_w(i*sp+1-k,t)
%                 nds(i*d+1-k)=1;
%             end
%         end
%     end
% 
% end








