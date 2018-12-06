clear variables; clc
% Script to demonstrate the issues around material objectivity described in
% the manuscript
%
% David Nolan
%
%% Material properties
props(1) = 0.2;   % C10 (MPa)
props(2) = 2;       % D1 (MPa^-1)
props(3) = 0.01;    % k1 (MPa)
props(4) = 2;     % k2 (-)
props(5) = 20;      % theta (degrees)
props(6) = -20;      % theta (degrees)

%% Compute the global and local deformation gradients
%
% Deformation gradient in global basis, FG
% N.B. This is what Abaqus passes to the UMAT when no local coordinate
% system is implemented.
FG1 = [1.10   0.10   0.00
      0.05   0.90   0.15
      0.20   0.00   1.20];      % Eq. (4)

% Compute a polar decomposition of the deformation gradient
[UG1, RG1]=polardecomp(FG1);

% Calculate the local deformation gradient
% FL is what Abaqus returns to the UMAT when you use *orientation
FL1 = RG1'*FG1*RG1;

% Get polar decomp of this local dgrad
[UL1, RL1]=polardecomp(FL1);

%% Calculate the stress in the global coordinate system
[sigGlob1] = MA_global(FG1,props);
pGlob1 = -trace(sigGlob1)/3;

%% Calculate the stress in the local coordinate system
[sigLoc1] = MA_local(FL1,props);
pLoc1 = -trace(sigLoc1)/3;

%% Calculate the co-rotational stress
[sigGlobU1] = MA_global(UG1,props);

%% Map the stress back from the local to the global coordinate system
stressTrans1 = RG1*sigLoc1*RG1';


%%
disp('########## Before Rotation ############')
disp('**** Stress Invariants ****')
disp('- Eigenvalues of global stress (Principal stresses)')
disp(eigs(sigGlob1))
disp('- Pressure stress of global stress')
disp(pGlob1)
disp('- Eigenvalues of local stress (Principal stresses)')
disp(eigs(sigLoc1))
disp('- Pressure stress of local stress')
disp(pLoc1)

disp('**** Stress Tensors ****')
disp('- Global stress based on global deformation gradient: FG')
disp(sigGlob1)
disp('- Local stress based on local deformation gradient: FL')
disp(sigLoc1)
disp('- Co-rotational stress based on right stretch tensor: UG')
disp(sigGlobU1)
disp('- Local stress transformed back to the global basis: R*sigLoc*R^T')
disp(stressTrans1)

%% Perform a rigid body rotation on the stretched configuration above

% Create rotation matrix Phi to rotate -pi/6 about the z axis
t=-pi/6;
Phi = [cos(t) -sin(t) 0
       sin(t)  cos(t) 0
       0       0      1];

% Rotate the stretch above by pi/6
FG2 = Phi*FG1;

% Now follow the same proceduce as above

% Compute a polar decomposition of the deformation gradient
[UG2, RG2]=polardecomp(FG2);

% Calculate the local deformation gradient
% FL is what Abaqus returns to the UMAT when you use *orientation
FL2 = RG2'*FG2*RG2;


%% Calculate the stress in the global coordinate system
[sigGlob2] = MA_global(FG2,props);
pGlob2 = -trace(sigGlob2)/3;

%% Calculate the stress in the local coordinate system
[sigLoc2] = MA_local(FL2,props);
pLoc2 = -trace(sigLoc1)/3;

%% Calculate the co-rotational stress
[sigGlobU2] = MA_global(UG2,props);
stressTrans2 = RG2*sigLoc2*RG2';

%%
disp('______________________________________________________________')
disp(' ')
disp('########## After Rotation ############')
disp('**** Stress Invariants ****')
disp('- Eigenvalues of global stress (Principal stresses)')
disp(eigs(sigGlob2))
disp('- Pressure stress of global stress')
disp(pGlob2)
disp('- Eigenvalues of local stress (Principal stresses)')
disp(eigs(sigLoc2))
disp('- Pressure stress of local stress')
disp(pLoc1)

disp('**** Stress Tensors ****')
disp('- Global stress based on global deformation gradient: FG')
disp(sigGlob2)
disp('- Local stress based on local deformation gradient: FL')
disp(sigLoc2)
disp('- Co-rotational stress based on right stretch tensor: UG')
disp(sigGlobU2)
disp('- Local stress transformed back to the global basis: R*sigLoc*R^T')
disp(stressTrans2)