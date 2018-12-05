clear variables; clc
% Script to demonstrate the issues around material objectivity described in
% the manuscript
%
% David Nolan
%
%%
% Material properties
props(1) = 0.2;   % C10 (MPa)
props(2) = 2;       % D1 (MPa^-1)
props(3) = 0.05;    % k1 (MPa)
props(4) = 2;     % k2 (-)
props(5) = 20;      % theta (degrees)
props(6) = -20;      % theta (degrees)
%%
% Deformation gradient in global basis, FG
% N.B. This is what Abaqus passes to the UMAT when no local coordinate
% system is implemented.
FG = [1.10   0.10   0.00
      0.05   0.90   0.15
      0.20   0.00   1.20];      % Eq. (4)

% Compute a polar decomposition of the deformation gradient
[UG, RG]=polardecomp(FG);

% Calculate the local deformation gradient
% FL is what Abaqus returns to the UMAT when you use *orientation
FL = RG'*FG*RG;

% Get polar decomp of this local dgrad
[UL, RL]=polardecomp(FL);

%% Calculate the stress in the global coordinate system
[sigGlob] = MA_global(FG,props);
pGlob = -trace(sigGlob)/3;


%% Calculate the stress in the local coordinate system
[sigLoc] = MA_local(FL,props);

pLoc = -trace(sigLoc)/3;

%% Calculate stress using
[sigGlobU] = MA_global(UG,props);
[sigLocU] = MA_local(UL,props);

stressTrans = RG*sigLoc*RG';


%%
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
disp('- Local stress transformed back to the global basis: R*sigLoc*R^T')
disp(stressTrans)
