% mat = rpy_ewald_matrix_mex(pos, rdi, L, tol, mode, ...)
%
% Compute mobility matrix using Ewald sum of RPY tensor.
%
% INPUT
%   pos     the array of particle coordinates
%   rdi     the array of particle radii
%   L       the simulation box size
%   tol     the requested tolerance of Ewald errors
%   mode    'full', 'real' or 'recip'
%   xi      the Ewald parameter (optional);
%           if not specified, the a value of xi will be
%           automatically chosen by ...
%
% OUTPUT
%   mat     the array of the mobility matrix
