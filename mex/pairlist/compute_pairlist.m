function [pi, pj] = compute_verlet(pos, rdi, L, cutoff, mode)
% [pi, pj] = compute_verlet(pos, rdi, L, cutoff, mode)
%   Compute list of all pairs of particles within a distance "cutoff"
%   of each other.  Periodic boundary conditions are used.
%   Output particle indices use base 0 and range from 0 to npos-1.
%   Note: algorithm will return "interaction" with itself.
%   Algorithm uses Verlet cell lists.  Cell size is chosen
%   based on cutoff and radius of largest particle.
%
% INPUT
%   pos     the array of particle coordinates (3*npos x 1)
%   rdi     the array of particle radii (npos x 1)
%   L       the simulation box size (scalar)
%   cutoff  the cutoff (scalar)
%   mode    'R' (absolute distance) or 'S' (normalized distance)
%
% OUTPUT
%   pi      the vector of particle indices (first of pair)
%   pj      the vector of particle indices (second of pair)
