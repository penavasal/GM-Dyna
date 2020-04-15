function main_plot(driver,varargin)


% Driver to employ plot tools developed within this code
% Input
%     driver:    - 'VTK' write *.vtk files to plot in paraview
%                - 'GRID' use griddata and surface to plot 

    VERSION='Code';
    RUTA_unix='~/REPOs/GM-Dyna/';
    RUTA_win ='\REPOs\GM-Dyna/';

    % PATH to files
    if ismac || isunix  % Code to run on Mac or Unix plaform 
        s=strcat(RUTA_unix,VERSION);
    elseif ispc         % Code to run on Windows platform
        s=strcat(getenv('HOMEPATH'),RUTA_win,VERSION);
    else
        disp('Platform not supported')
        stop
    end 

    addpath(path,s);
    folder  = pwd;
    PLOT(driver,folder,varargin{:});
    


    
