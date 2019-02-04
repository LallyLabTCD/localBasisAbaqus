clear variables

%% Continuum element where E_i = G_i
% Basic deformation
F = [1 0 0
     0 1.1 0
     0 0 1/1.1];

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
 
%% Continuum element where E_i \= G_i
% Rotation matrix to rotate basis in reference configuration by 20 degrees
t=20;
Q = [cosd(t) -sind(t) 0
     sind(t)  cosd(t) 0
       0         0    1];

% "Abaqus Local" deformation gradient where the reference basis is not
% aligned with the global basis
Fal20 = Q'*Fal*Q;

%% Shell element where E_i \= G_i
% Rotation matrix to rotate basis
t=20;
Q = [cosd(t) -sind(t) 0
     sind(t)  cosd(t) 0
       0         0    1];

% "Classical" local deformation gradient where the reference basis is not
% aligned with the global basis
Fd20 = Q'*Fd*Q;

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
disp('***** Continuum Element: E_i neq G_i *****')
disp(' ')
disp('Fal20 = ')
disp(Fal20)
disp(' ')
disp('***** Structural Element: E_i neq G_i *****')
disp(' ')
disp('Fd20 = ')
disp(Fd20)
disp(' ')
   
   