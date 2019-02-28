clear variables
% Script to demonstrate the difference between structural elements (e.g. 
% shells and membranes) and continuum elements (here, plane stress).
%
% David Nolan

%% Continuum element where E_i = G_i
% Basic deformation
F = [1 0   0
     0 1.1 0
     0 0   1/1.1];

% Rotate by 45 degrees
t=45;
T = [cosd(t)   sind(t) 0
     -sind(t)  cosd(t) 0
       0         0    1];

% Rotate by 45 degrees. DG in the global basis
F = T*F;
disp(F)

% Determine rotational part of the deformation.
[U,R]=polardecomp(F);

% "Abaqus Local" deformation gradient
Fal = R'*F*R;

%% Shell element where E_i = G_i

% Polar decomposition to determine R  
[U,R]=polardecomp(F);

% "Classical" local deformation gradient
Fd = R'*F;
 

%% Display Results
clc
disp('***** Continuum Element: Global Basis *****')
disp(' ')
disp('F = ')
disp(F)
disp(' ')
disp('***** Continuum Element: Abaqus Local Basis *****')
disp(' ')
disp('Fal = ')
disp(Fal)
disp(' ')
disp('***** Structural Element *****')
disp(' ')
disp('Fd = ')
disp(Fd)
disp(' ')
   
   