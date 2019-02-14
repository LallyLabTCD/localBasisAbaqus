clear variables

%% Continuum element where E_i = G_i
% Basic deformation
F = [1 0 0
     0 2 0
     0 0 0.5];

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
disp(' +++ N.B. this is not the same as the F returned by WE1_PStress_Ori_largeDef.inp ')
disp(' ')
disp('***** Structural Element *****')
disp(' ')
disp('Fd = ')
disp(Fd)
disp(' ')
   