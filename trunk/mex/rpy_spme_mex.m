% v = rpy_spme_mex(pos, rdi, L, f, mode, xi, rmax, dim, porder) 
%
% Compute product of mobility matrix with force vector using SPME.
%
% INPUT
%   pos       particle coordinates (3 x npos)
%   rdi       particle radii (1 x npos)
%   L         simulation box width (cubical L x L x L box)
%   f         forces (3npos x nrhs)
%   mode      'full', 'real' or 'recip'
%   xi        Ewald parameter
%   rmax      real-space cutoff
%   dim       FFT grid dimension
%   porder    interpolation order
%
% OUTPUT
%   v         velocities (3npos x nrhs)
