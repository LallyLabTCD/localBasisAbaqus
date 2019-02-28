clear variables; 
% Script to demonstrate the issues around material objectivity described in
% the manuscript. This is an isotropic material model example where E \neq G
%
% David Nolan
%
% Material properties
props(1) = 0.2;   % C10 (MPa)
props(2) = 2;     % D1 (MPa^-1)
%% Kinematics
% Deformation gradient in global basis, FG
% N.B. This is what Abaqus passes to the UMAT when no local coordinate
% system is implemented.
FG = [1.10   0.10   0.00
      0.05   0.90   0.15
      0.20   0.00   1.20];      % Eq. (4);

% Rotation matrix used to change the basis in the reference
% configuration.
t=0.45; % 20 degrees in radians
Q = [cosd(t) -sind(t) 0
     sind(t)  cosd(t) 0
       0     0      1];         % Eq. (7)

% Change the basis of the global deformation gradient
FG20 = Q'*FG*Q;                 % Eq. (8)

% Polar decomposition of the global deformation gradient in the new basis
[UG20, RG20]=polardecomp(FG20);

% Map the global deformation gradient to the local deformation gradient
FL20 = RG20'*FG20*RG20;         % Eq. (9)
[UL20, RL20]=polardecomp(FL20);

%% Calculate the stress
% Calculate the stress in the global coordinate system. Eq. (11)
[sigGlob] = NeoHooke(FG,props);
% Calculate pressure stress
pGlob = -trace(sigGlob)/3;

% Calculate the stress in the local coordinate system. Term in curly
% brackets in Eq. (14)
[sigLoc] = NeoHooke(FL20,props);
% Calculate pressure stress
pLoc = -trace(sigLoc)/3;

% Calculate the co-rotational stress. Term in curly brackets in Eq. (15)
[sigGlobU] = NeoHooke(UG20,props);

% Change the basis of the local stress back to the global basis. Eq (19)
stressTrans = Q*RG20*sigLoc*RG20'*Q';

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