% Code to demonstrate the calculation of the local deformation gradient
% from the global deformation gradient  

clear variables
% Deformation gradient
FG = [1.10   0.10   0.00
      0.05   0.90   0.15
      0.20   0.00   1.20];      % Eq. (22)

% Global basis vectors
G1 = [1 0 0]';
G2 = [0 1 0]';
G3 = [0 0 1]';

% The local basis vectors in the reference configuration. Note that in this 
% case G_i = E_i
E1 = G1;
E2 = G2;
E3 = G3;

% Store basis vectors in this array
E=[E1 E2 E3];


%% Global deformation gradient
% e_i = E_i = G_i 
%
% This section only serves to show how the index notation equations are
% implmented using for loops. F will equal FG here.
F=zeros(3);
for i=1:3
    for j=1:3
        F = F + FG(i,j)*E(:,i)*E(:,j)';
    end
end


%% Determine Co-rotation R through polar decomposition
% R - rotation matrix
% U - right stretch tensor
[U, R]=polardecomp(F);


%% Calculate the local basis vectors in the current configuration.
% Co-rotational basis system

eL1 = R*E1;
eL2 = R*E2;
eL3 = R*E3;
% Store basis vectors in this array
e=[eL1 eL2 eL3];


%% Formulate a Q matrix
% Merely a demonstration of how this is calculated from basis vectors
% We will find that Qij = R;

Qij=zeros(3);   % Dummy matrix
for i=1:3
    for j=1:3
        Qij(i,j) = dot(E(:,i), e(:,j));
    end
end


%% Calculation of the local deformation gradient (Eq. 65)

Fd = zeros(3);  % Dummy matrix

% Matrix multiplication using indices
for i=1:3
    for j=1:3
        for k=1:3
            Fd(j,k) = Fd(j,k) + FG(i,j)*Qij(i,k);
        end
    end
end

% Alternative expression
Fd_alt = F'*Qij;


%% Use the local deformation gradient to calculate the global 

Fnew=zeros(3);  % Dummy matrix
for k=1:3
    for j=1:3
        Fnew = Fnew + Fd(j,k)*e(:,k)*E(:,j)';
    end
end
