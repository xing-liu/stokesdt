function [pi, pj] = compute_verlet(pos, rdi, L, cutoff, mode)
% [pi, pj] = compute_verlet(pos, rdi, L, cutoff, mode)
%   Compute short-range interactions between all the particles.
%   The particle pairs within the cutoff distance are stored in [pi, pj].
% INPUT
%   pos     the array of particle coordinates (3*npos x 1)
%   rdi     the array of particle radii (npos x 1)
%   L       the simulation box size
%   cutoff  the cutoff
%   mode    'R' (absolute distance) or 'S' (normalized distance)
%
% OUTPUT
%   pi      the vector of particle indices (first of pair)
%   pj      the vector of particle indices (second of pair)
