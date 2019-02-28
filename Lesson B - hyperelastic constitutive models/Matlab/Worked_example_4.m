clear variables; clc
% Script to demonstrate the issues around material objectivity described in
% the manuscript. This is an example using an isotropic material model.
%
% David Nolan
%
% Material properties for the neo-Hookean constitutive model
props(1) = 0.2;   % C10 (MPa)
props(2) = 2;     % D1 (MPa^-1)
%% Kinematics
% Deformation gradient in global basis, FG
% N.B. This is what Abaqus passes to the UMAT when no local coordinate
% system is implemented.
FG = [1.10   0.10   0.00
      0.05   0.90   0.15
      0.20   0.00   1.20];      % Eq. (4)

% Compute a polar decomposition of the deformation gradient
[UG, RG]=polardecomp(FG);

% Calculate the local deformation gradient
% N.B. This is what Abaqus passes to the UMAT when you use *orientation and
% E_i = G_i
Fal = RG'*FG*RG;

% Get polar decomp of this local deformation gradient
[UL, RL]=polardecomp(Fal);

%% Calculate the stress
% Calculate the stress in the global coordinate system. Eq. (11)
[sigGlob] = NeoHooke(FG,props);
% Calculate pressure stress
pGlob = -trace(sigGlob)/3;

% Calculate the stress in the local coordinate system. Term in curly
% brackets in Eq. (14)
[sigLoc] = NeoHooke(Fal,props);
% Calculate pressure stress
pLoc = -trace(sigLoc)/3;

% Calculate the co-rotational stress. Term in curly brackets in Eq. (15)
[sigGlobU] = NeoHooke(UG,props);

% Change the basis of the local stress back to the global basis. Eq (14)
stressTrans = RG*sigLoc*RG';

%% Output
disp('**** Stress Invariants ****')
disp('- Eigenvalues of global stress (Principal stresses)')
disp(eigs(sigGlob))
disp('- Pressure stress of global stress')
disp(pGlob)
disp('- Eigenvalues of local stress (Principal stresses)')
disp(eigs(sigLoc))
disp('- Pressure stress of local stress')
disp(pLoc)

disp('**** Stress Tensors ****')
disp('- Global stress based on global deformation gradient: FG')
disp(sigGlob)
disp('- Local stress based on local deformation gradient: FL')
disp(sigLoc)
disp('- Co-rotational stress based on right stretch tensor: UG')
disp(sigGlobU)
disp('- Local stress transformed back to the global basis')
disp(stressTrans)