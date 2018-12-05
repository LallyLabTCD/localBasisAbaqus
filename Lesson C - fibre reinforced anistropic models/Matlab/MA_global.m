function [cstress, kirchIso, kirchAniso, kirch] = MA_global(F, props)
%MA_GLOBAL Calculate stress defined by MA hyperelastic constitutive model
%   MA_GLOBAL(F, props) returns the stress for a give deformation gradient
%   F and set of material properties PROPS.

%% Preliminaries
% Import material properties
C10 = props(1);
D1 = props(2);
k1 = props(3);
k2 = props(4);
thetad1 = props(5);
thetad2 = props(6);
%
% Convert degrees to radians
thetar1 = thetad1*pi/180;
thetar2 = thetad2*pi/180;

%% Fibre kinematics
% Create fibre orientation vectors in the undeformed configuration
a01 = [cos(thetar1) sin(thetar1)  0]';
a02 = [cos(thetar2) sin(thetar2) 0]';

% Fibre orientation vector in the deformed configuration
a1 = F*a01; 
a2 = F*a02;

% Collect the fibre vectors in an array
Amat = [a1 a2];


%% Calculate the isotropic stress
% Invariants (it doesn't matter if the invarints are calculated using
% the left or right deformation tensor (C or B) )
B = F*F';
J = det(F);
J23 = J^(2/3);
I1 = trace(B);

% The volumetric part of the isotropic stress
kirchIso = 2*(J-1)*J/D1*eye(3,3);

% The isochoric part of the isotropic stress
kirchIso = kirchIso + 2*C10*B'/J23 - 2*C10*I1/3/J23*eye(3,3);


%% Calcaulate the anisotropic portion
% Initialize the anisotropic stress
kirchAniso = zeros(3,3);

% Loop over the number of fibres in the model
for i = 1:2
    
    % Extract the i'th fibre vector
    ai = Amat(:,i);
    
    % Compute the dyadic product of the deformed orientation vectors
    axa = ai*ai';
    
    % Compute the anisotropic deformation invariant
    I4 = trace(axa);
    
    % The fibre will only contribute if it is in tension. i.e. if its
    % stretch is greater than 1 (N.B. I4 is the square of fibre stretch)
    if (I4 > 1)
        
        % The anisotropic term is broken into 3 parts to avoid mistakes
        % with brackets & multiplication
        p1 = 2*k1*(I4-1);

        p2 = exp(k2*((I4-1)^2));

        p3 = axa;
        
        % Add the contribution of the i'th fibre to the anisotropic portion
        % of the stress
        kirchAniso = kirchAniso + p1*p2*p3;
    
    end
    
end

% To calculate total stress, sum the isotropic and anisotropic parts
kirch = kirchIso + kirchAniso;

% Convert the Kirchhoff stress to Cauchy stress
cstress = kirch/J;