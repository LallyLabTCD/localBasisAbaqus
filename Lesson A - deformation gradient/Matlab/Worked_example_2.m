clear variables; clc
% Script to demonstrate the issues around material objectivity described in
% the manuscript
%
% David Nolan
%
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
% This is the same transform that Abaqus uses when you use *orientation
Fal = RG'*FG*RG;                 % Eq. (3) & (6)
[UL, RL]=polardecomp(Fal);       % Eq. (1)

% Note the the eigenvalues of the right stretch tensor U are the same. We
% expect to see this invariance.
[lambdaGvec, lambdaGval] = eigs(UG);
[lambdaLvec, lambdaLval] = eigs(UL);

%% Basis vectors
% Global basis vectors e0G
G1 = [1 0 0]';
G2 = [0 1 0]';
G3 = [0 0 1]';

% The local basis vectors in the reference configuration. Note that in this 
% case G_i = E_i
E1 = G1;
E2 = G2;
E3 = G3;

% Calculate the local basis vectors in the current configuration.
eL1 = RG*E1;
eL2 = RG*E2;
eL3 = RG*E3;


%% Plotting
figure(1);clf;
title('Local Basis System')
% Plot the deformed unit cube
dfgrdplot(FG,'w')

% Plot and label the global basis vectors
quiver3(0,0,0,G1(1),G1(2),G1(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,G2(1),G2(2),G2(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,G3(1),G3(2),G3(3),'k','LineWidth',2,'MaxHeadSize',0.4)

G = [G1, G2 G3];
strings = {'$$\textbf{G}_{1}$$';'$$\textbf{G}_{2}$$';'$$\textbf{G}_{3}$$'}; 
text(G(1,:),G(2,:),G(3,:),strings,'Interpreter','latex','FontSize',14)

% Plot and label the local basis vectors in the current configuration
quiver3(0.5,0.5,0.5,eL1(1),eL1(2),eL1(3),'r','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0.5,0.5,0.5,eL2(1),eL2(2),eL2(3),'g','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0.5,0.5,0.5,eL3(1),eL3(2),eL3(3),'b','LineWidth',2,'MaxHeadSize',0.4)

e0L = [eL1, eL2 eL3] + 0.5;
strings = {'$$\textbf{e}_{1}$$';'$$\textbf{e}_{2}$$';'$$\textbf{e}_{3}$$'}; 
text(e0L(1,:),e0L(2,:),e0L(3,:),strings,'Interpreter','latex','FontSize',14)

% Set axis limits
xlim([-0.2 1.5])
ylim([-0.2 1.5])
zlim([-0.2 1.5])
 
% Label the basis vectors
e0L = [eL1, eL2 eL3] + 0.5;
strings = {'$$\textbf{e}_{1}$$';'$$\textbf{e}_{2}$$';'$$\textbf{e}_{3}$$'}; 
text(e0L(1,:),e0L(2,:),e0L(3,:),strings,'Interpreter','latex','FontSize',14)

% Set the camera view
view(125,25)