function dfgrdplot(F, colour)
% DFGRDPLOT plot the deformation gradient
%   
%   DFGRDPLOT(F, colour) - This function will plot the deformation of a cube 
%   that will occur when a prescribed deformation gradient F is applied to 
%   it. F must be a 3x3 matrix, and is the identity when there is no 
%   deformation. The cube is semi-transparent and its colour may be changed
%   using string variables such as 'r' for red etc.

% Matrix containing the coordinates for each vertex of the cube
orig_coords = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
orig_coords = orig_coords';
[l, b] = size(orig_coords);     % Dimension of the array
x = zeros(l,b);                 % preallocation of deformed coord

% Subject each vertex to the deformation prescribed by the deformation
% gradient
    for j = 1:b        
        X = orig_coords(:,j);   % Original coordinate
        x(:,j) = F*X;           % Deformed coordinate
    end
    
    
% Connectivity of the faces to the vertices
f = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];


% Patch command will plot the cube
hold on
patch('Vertices',x','Faces',f, 'FaceColor',colour,'FaceAlpha',0.5);

% Setup for the figure
view(130., 30);
grid on
axis square  
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontSize',16)
% axis([-0.5 1.5 0 2 0 2]);