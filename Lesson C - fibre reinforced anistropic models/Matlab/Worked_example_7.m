clear variables; clc
% Script to demonstrate the issues around material objectivity described in
% the manuscript
%
% David Nolan
%
%% Global kinematics

% Deformation gradient in global basis, FG
% N.B. This is what Abaqus passes to the UMAT when no local coordinate
% system is implemented.
FG = [1.10   0.10   0.00
      0.05   0.90   0.15
      0.20   0.00   1.20];      % Eq. (4);

% Fibre vector in the reference configuration in the global coordinate
% frame
A_G = [cosd(30) sind(30) 0]';

% Fibre vector in the current configuration in the global coordinate
% frame
a_G = FG*A_G;

% The structural tensor in the global coordinate frame
aoaG = a_G*a_G';

%% Local kinematics

% Rotation matrix used to change the basis in the reference
% configuration.
t=20; % 20 degrees in degrees
Q = [cosd(t) -sind(t) 0
     sind(t)  cosd(t) 0
       0         0    1];         % Eq. (7)

% Change the basis of the global deformation gradient
FG20 = Q'*FG*Q;                   % Eq. (8)
% Change the basis of the global structural tensor
aoaG20 = Q'*aoaG*Q;

% Polar decomposition of the global deformation gradient in the new basis
[UG20, RG20]=polardecomp(FG20);

% Map the global deformation gradient to the local deformation gradient
Fal20 = RG20'*FG20*RG20;         % Eq. (9)

% Define the fibre vector in the reference configuration with respect to
% the local coordinate basis vectors E_i.
%
% As the global fibre angle is 30deg with respect to G_1 and the local
% E_1 vector is at 20deg with respect to G_1, we choose a fibre angle of
% 30-20 = 10 degrees, to ensure equivalence of the local and global cases.
A_E = [cosd(10) sind(10) 0]';

% Convert the Abaqus Local deformation gradient to U
U20 = Fal20*RG20';

% Deformed fibre vector in the current configuration in the local basis
% system e_i
a_e = U20*A_E;

% Calculate the structural tensor using the fibre vector in the current
% configuration 
aoaL20 = Fal20*RG20'*A_E*A_E'*RG20*Fal20';

% Alternatively: aoaL20 = a_e*a_e'


% Calculate the anisotropic invariant I4 and verify that they are identical
% in all cases.
IfG = trace(aoaG);
IfG20 = trace(aoaG20);
IfL20 = trace(aoaL20);

%% Using the right stretch tensor
aUG_e = UG20*A_E;

aoaUG = aUG_e*aUG_e';

IfGL = trace(aoaUG);

%% Map the local structural tensor back to the global basis
aoaTrans = Q*RG20*aoaUG*RG20'*Q';

%% Output results
disp('**** Fibre Invariants ****')
disp('- I_f based in global basis')
disp(IfG20)
disp('- I_f based in local basis')
disp(IfL20)
disp('- I_f based on right stretch in global basis')
disp(IfGL)

disp('**** Structural tensor ****')
disp('- aoa based in global basis')
disp(aoaG)
disp('- aoa based in local basis')
disp(aoaL20)
disp('- aoa based on right stretch in global basis')
disp(aoaUG)
disp('- local aoa mapped back to global')
disp(aoaTrans)

%% Basis vectors
% Global basis vectors
G1 = [1 0 0]';
G2 = [0 1 0]';
G3 = [0 0 1]';

% The basis vectors in the reference configuration.
% Global basis
E1G = G1;
E2G = G2;
E3G = G3;

% Local basis
E1L = Q*G1;
E2L = Q*G2;
E3L = Q*G3;

% Calculate the basis vectors in the current configuration.
% Global basis
e1G = E1G;
e2G = E2G;
e3G = E3G;

% Local basis
e1L = RG20*E1L;
e2L = RG20*E2L;
e3L = RG20*E3L;

%% Plot the global scenario
figure(1);clf;
% Plot the deformed unit cube
dfgrdplot(FG,'w')

%---------

% Plot and label the global basis vectors
quiver3(0,0,0,G1(1),G1(2),G1(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,G2(1),G2(2),G2(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,G3(1),G3(2),G3(3),'k','LineWidth',2,'MaxHeadSize',0.4)
G = [G1, G2 G3];
strings = {'$$\textbf{G}_{1}$$';'$$\textbf{G}_{2}$$';'$$\textbf{G}_{3}$$'}; 
text(G(1,:),G(2,:),G(3,:),strings,'Interpreter','latex','FontSize',14)


% Plot and label the basis vectors in the current configuration.
% These are offset by 0.5 for display purposes only
quiver3(0.5,0.5,0.5,e1G(1),e1G(2),e1G(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0.5,0.5,0.5,e2G(1),e2G(2),e2G(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0.5,0.5,0.5,e3G(1),e3G(2),e3G(3),'k','LineWidth',2,'MaxHeadSize',0.4)
eG = [e1G, e2G e3G] + 0.5;
strings = {'$$\textbf{e}_{1}$$';'$$\textbf{e}_{2}$$';'$$\textbf{e}_{3}$$'}; 
text(eG(1,:),eG(2,:),eG(3,:),strings,'Interpreter','latex','FontSize',14)

%---------

% Plot and label the fibre vector in the undeformed configuration
quiver3(0,0,0,A_G(1),A_G(2),A_G(3),'b','LineWidth',2,'MaxHeadSize',0.4)
text(A_G(1),A_G(2),A_G(3),'$$\textbf{A}_{G_i}$$','Interpreter','latex','FontSize',16,'Color','b')

% Plot the fibre vector in the deformed configuration
quiver3(0.5,0.5,0.5,a_G(1),a_G(2),a_G(3),'b','LineWidth',2,'MaxHeadSize',0.4)
text(a_G(1)+0.5,a_G(2)+0.5,a_G(3)+0.5, '$$\textbf{a}_{G_i}$$','Interpreter','latex','FontSize',16,'Color','b')

%---------

% Set axis limits
xlim([-0.4 1.5])
ylim([-0.4 1.5])
zlim([-0.4 1.5])
% Set the camera view
view(0,90)

%% Plot the local scenario
figure(2);clf;
% Plot the deformed unit cube
dfgrdplot(FG,'w')

%---------

% Plot and label the global basis vectors
quiver3(0,0,0,G1(1),G1(2),G1(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,G2(1),G2(2),G2(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,G3(1),G3(2),G3(3),'k','LineWidth',2,'MaxHeadSize',0.4)
G = [G1, G2 G3];
strings = {'$$\textbf{G}_{1}$$';'$$\textbf{G}_{2}$$';'$$\textbf{G}_{3}$$'}; 
text(G(1,:),G(2,:),G(3,:),strings,'Interpreter','latex','FontSize',14)

% Plot and label the local basis vectors in the undeformed config.
quiver3(0,0,0,E1L(1),E1L(2),E1L(3),'Color',[0.87 0.49 0],'LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,E2L(1),E2L(2),E2L(3),'Color',[0.87 0.49 0],'LineWidth',2,'MaxHeadSize',0.4)
quiver3(0,0,0,E3L(1),E3L(2),E3L(3),'Color',[0.87 0.49 0],'LineWidth',2,'MaxHeadSize',0.4)
E = [E1L, E2L E3L];
strings = {'$$\textbf{E}_{1}$$';'$$\textbf{E}_{2}$$';'$$\textbf{E}_{3}$$'}; 
text(E(1,:),E(2,:),E(3,:),strings,'Interpreter','latex','FontSize',14)

% Plot and label the local basis vectors in the current configuration
quiver3(0.5,0.5,0.5,e1L(1),e1L(2),e1L(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0.5,0.5,0.5,e2L(1),e2L(2),e2L(3),'k','LineWidth',2,'MaxHeadSize',0.4)
quiver3(0.5,0.5,0.5,e3L(1),e3L(2),e3L(3),'k','LineWidth',2,'MaxHeadSize',0.4)
e0L = [e1L, e2L e3L] + 0.5;
strings = {'$$\textbf{e}_{1}$$';'$$\textbf{e}_{2}$$';'$$\textbf{e}_{3}$$'}; 
text(e0L(1,:),e0L(2,:),e0L(3,:),strings,'Interpreter','latex','FontSize',14)


% Plot the fibre vector in the undeformed configuration.
% N.B. because we are plotting with respect to a global coordinate system
% we use AG and aG here even though they represent the local basis vectors.
quiver3(0,0,0,A_G(1),A_G(2),A_G(3),'b','LineWidth',2,'MaxHeadSize',0.4)
text(A_G(1),A_G(2),A_G(3),'$$\textbf{A}_{E_i}$$','Interpreter','latex','FontSize',16,'Color','b')
% Plot the fibre vector in the deformed configuration
quiver3(0.5,0.5,0.5,a_G(1),a_G(2),a_G(3),'b','LineWidth',2,'MaxHeadSize',0.4)
text(a_G(1)+0.5,a_G(2)+0.5,a_G(3)+0.5, '$$\textbf{a}_{e_i}$$','Interpreter','latex','FontSize',16,'Color','b')

%---------

% Set axis limits
xlim([-0.4 1.5])
ylim([-0.4 1.5])
zlim([-0.4 1.5])
% Set the camera view
view(0,90)