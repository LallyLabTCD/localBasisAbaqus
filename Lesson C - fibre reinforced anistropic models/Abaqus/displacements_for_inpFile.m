%  F1 = [1.1 0.05 0
%     0.0 1.05 0
%     0 0 1];
% 
F1 = [1.1   0.1 0
      0.05   0.9   0.15
      0.2  0   1.2];

orig_coords = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
orig_coords = orig_coords';
[l, b] = size(orig_coords);     % Dimension of the array
x1 = zeros(l,b);                 % preallocation of deformed coord

% Subject each vertex to the deformation prescribed by the deformation
% gradient
    for j = 1:b        
        X = orig_coords(:,j);   % Original coordinate
        x1(:,j) = F1*X;           % Deformed coordinate
    end
U1 = x1-orig_coords;
    
%%
t=-pi/6;
Rc = [cos(t) -sin(t) 0
       sin(t) cos(t) 0
       0 0 1];
   
F2 = Rc*F1;
x2 = zeros(l,b);                 % preallocation of deformed coord
for j = 1:b        
    X = orig_coords(:,j);   % Original coordinate
    x2(:,j) = F2*X;           % Deformed coordinate
end

U2 = x2-orig_coords;
 
% figure(2);clf;hold all
% plot3(orig_coords(1,:),orig_coords(2,:),orig_coords(3,:),'bo')
% plot3(x1(1,:),x1(2,:),x1(3,:),'ro')
% plot3(x2(1,:),x2(2,:),x2(3,:),'yo')
% axis equal

 