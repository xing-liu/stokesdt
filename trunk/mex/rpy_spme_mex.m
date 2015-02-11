% v = rpy_spme_mex(pos, rdi, L, f, mode, xi, rmax, dim, porder) 
%
% Compute product of mobility matrix with force vector using SPME.
%
% INPUT
%   pos       the vector of particle coordinates
%   rdi       the vector of particle radii
%   L         the simulation box size
%   f         the block of vectors of forces
%   mode      'full', 'real' or 'recip'
%   xi        the Ewald parameter
%   rmax      the real-space cutoff
%   dim       the dimension of FFT grid
%   porder    the interpolation order
% OUTPUT
%   v         the block of vectors of velocities
