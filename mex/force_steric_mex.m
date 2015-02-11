% f = force_steric_mex(pos, rdi, L, r0, k0)
%
% Compute steric forces given particle positions.
%
% INPUT
%   pos       particle coordinates (3 x npos)
%   rdi       particle radii (1 x npos)
%   L         simulation box width (cubical L x L x L box)
%   r0        cutoff distance (often chosen as 2)
%   k0        force constant (often chosen as 125)
%
% OUTPUT
%   f         forces (3 x npos)
