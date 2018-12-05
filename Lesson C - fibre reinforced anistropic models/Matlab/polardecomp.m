function [U, R] = polardecomp(F)
%POLARDECOMP Compute the polar decomposition of a deformation gradient.
%   [U, R] = POLARDECOMP(F) takes the deformation gradient F and decomposes
%   it in into a rigid rotation part R and the right stretch tensor U. 
%   The multiplicative decomposition is defined as F = RU

%%
C = F'*F; % right Cauchy Green Tensor

% Right stretch tensor
[Cvec, Cval] = eigs(C);
U = Cvec*sqrt(Cval)/Cvec;

% Rotation in undef config
R = F/U;