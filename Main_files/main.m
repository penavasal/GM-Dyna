function main(varargin)

%clear all
clearvars -except varargin
format long
restoredefaultpath
warning('off', 'MATLAB:nearlySingularMatrix')

% Main file
% Version Dic 2019

VERSION='Code';
RUTA_unix='~/REPOs/GM-Dyna/';
RUTA_win ='\REPOs\GM-Dyna/';
 
% PATH to files
if ismac || isunix  % Code to run on Mac or Unix plaform 
    s=strcat(RUTA_unix,VERSION);
    s1=strcat(s,'/Model');
elseif ispc         % Code to run on Windows platform
    s=strcat(getenv('HOMEPATH'),RUTA_win,VERSION);
    s1=strcat(s,'\Model');
else
    disp('Platform not supported')
    stop
end 

addpath(path,s);
addpath(path,s1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE PARAMETERS, GEOMETRY, MATERIAL and BOUNDARY CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0
    str='problem.txt';
else
    str=strcat(varargin{1},'.txt');
end
[MAT_POINT]=read_problem(str);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SOLVER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOLVER(MAT_POINT);



