clear variables; clc
% Script to demonstrate the issues around material objectivity described in
% the manuscript
%
% David Nolan
%
%%

% Deformation gradient in global basis, FG
% N.B. This is what Abaqus passes to the UMAT when no local coordinate
% system is implemented.
FG = [1.10   0.10   0.00
      0.05   0.90   0.15
      0.20   0.00   1.20];      % Eq. (4);

% Rotation matrix used to change the basis in the reference
% configuration.
t=pi/9; % 20 degrees in radians
Q = [cos(t) -sin(t) 0
     sin(t)  cos(t) 0
       0     0      1];         % Eq. (7)

% Change the basis of the global deformation gradient
FG20 = Q'*FG*Q;                 % Eq. (8)

% Polar decomposition of the global deformation gradient in the new basis
[UG20, RG20]=polardecomp(FG20);

% Map the global deformation gradient to the local deformation gradient
Fal20 = RG20'*FG20*RG20;         % Eq. (9)
[Ual20, RL20]=polardecomp(Fal20);

% Note the the eigenvalues of the right stretch tensor U are the same. We
% expect to see this invariance.
[lambdaGvec, lambdaGval] = eigs(UG20);
[lambdaLvec, lambdaLval] = eigs(Ual20);

%% Basis vectors
% Global basis vectors e0G
G1 = [1 0 0]';
G2 = [0 1 0]';
G3 = [0 0 1]';

% The basis vectors in the reference configuration.
E1 = Q*G1;
E2 = Q*G2;
E3 = Q*G3;

% Calculate the basis vectors in the current configuration.
eL1 = RG20*E1;
eL2 = RG20*E2;
eL3 = RG20*E3;


%% Plotting
figure(2);clf;
title('Local Basis System when E_i \neq G_i')
% Plot the deformed unit cube
dfgrdplot(FG,'w')

% Plot and label the global basis vectors
quiver3(0,0,0,G1(1),G1(2),G1(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,G2(1),G2(2),G2(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,G3(1),G3(2),G3(3),'k','LineWidth',2,'MaxHeadSize',0.4)

G = [G1, G2 G3];
strings = {'$$\textbf{G}_{1}$$';'$$\textbf{G}_{2}$$';'$$\textbf{G}_{3}$$'}; 
text(G(1,:),G(2,:),G(3,:),strings,'Interpreter','latex','FontSize',14)

% Plot and label the local basis vectors in the undefromed config.
quiver3(0,0,0,E1(1),E1(2),E1(3),'Color',[0.87 0.49 0],'LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,E2(1),E2(2),E2(3),'Color',[0.87 0.49 0],'LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,E3(1),E3(2),E3(3),'Color',[0.87 0.49 0],'LineWidth',2,'MaxHeadSize',0.4)

E = [E1, E2 E3];
strings = {'$$\textbf{E}_{1}$$';'$$\textbf{E}_{2}$$';'$$\textbf{E}_{3}$$'}; 
text(E(1,:),E(2,:),E(3,:),strings,'Interpreter','latex','FontSize',14)

% Plot and label the local basis vectors
quiver3(0.5,0.5,0.5,eL1(1),eL1(2),eL1(3),'r','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0.5,0.5,0.5,eL2(1),eL2(2),eL2(3),'g','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0.5,0.5,0.5,eL3(1),eL3(2),eL3(3),'b','LineWidth',2,'MaxHeadSize',0.4)

e0L = [eL1, eL2 eL3] + 0.5;
strings = {'$$\textbf{e}_{1}$$';'$$\textbf{e}_{2}$$';'$$\textbf{e}_{3}$$'}; 
text(e0L(1,:),e0L(2,:),e0L(3,:),strings,'Interpreter','latex','FontSize',14)

% Set axis limits
xlim([-0.4 1.5])
ylim([-0.4 1.5])
zlim([-0.4 1.5])
 
% Label the basis vectors
e0L = [eL1, eL2 eL3] + 0.5;
strings = {'$$\textbf{e}_{1}$$';'$$\textbf{e}_{2}$$';'$$\textbf{e}_{3}$$'}; 
text(e0L(1,:),e0L(2,:),e0L(3,:),strings,'Interpreter','latex','FontSize',14)

% Set the camera view
view(125,25)