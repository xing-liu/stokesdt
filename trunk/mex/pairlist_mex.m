% [ia, ja] = pairlist_mex(pos, rdi, L, cutoff, mode)
%
%   Compute list of all pairs of particles within a distance "cutoff"
%   of each other.  Periodic boundary conditions are used.
%   Output is "compressed sparse row" format giving a sparse matrix of pairs.
%   Note: algorithm will return "interaction" with itself.
%   Algorithm uses Verlet cell lists.  Cell size is chosen
%   based on cutoff and radius of largest particle.
%
% INPUT
%   pos     the array of particle coordinates (3 x npos)
%   rdi     the array of particle radii (1 x npos)
%   L       simulation box width (cubical L x L x L box)
%   cutoff  scalar cutoff distance
%   mode    'R' (absolute distance) or 'S' (normalized distance)
%
% OUTPUT
%   ia      row pointers, 0-based indexing
%   ja      column indices, 0-based indexing
