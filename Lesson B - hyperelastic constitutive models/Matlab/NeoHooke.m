function [cstress] = NeoHooke(F, props)
%NEOHOOKE Calculate stress defined by MA hyperelastic constitutive model
%   NEOHOOKE(F, props) returns the stress for a give deformation gradient
%   F and set of material properties PROPS.

%% Preliminaries
% Import material properties
C10 = props(1);
D1 = props(2);

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

cstress = kirchIso/J;