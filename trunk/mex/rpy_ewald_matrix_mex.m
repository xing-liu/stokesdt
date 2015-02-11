% mat = rpy_ewald_matrix_mex(pos, rdi, L, tol, mode, xi)
%
% Compute mobility matrix using Ewald sum of RPY tensor.
%
% INPUT
%   pos       particle coordinates (3 x npos)
%   rdi       particle radii (1 x npos)
%   L         simulation box width (cubical L x L x L box)
%   tol       Ewald error tolerance
%   mode      'full', 'real' or 'recip'
%   xi        Ewald parameter (optional) with default 10^(1/6)*pi/L)
%        
% OUTPUT
%   mat       mobility matrix (3npos x 3npos)
