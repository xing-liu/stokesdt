function [pi, pj] = compute_verlet(pos, rdi, L, cutoff, mode)
% [pi, pj] = compute_verlet(pos, rdi, L, cutoff, mode)
%   Compute short-range interactions between all the particles.
%   The particle pairs within the cutoff distance are stored in [pi, pj].
% INPUT
%   pos     the array of particle coordinates
%   rdi     the array of particle radii
%   L       the simulation box size
%   cutoff  the cutoff
%   mode    'R' (absulote distance) or 'S' (normalized distance)
%
% OUTPUT
%   pi      the vector of particle identities of pi
%   pj      the vector of particle identities of pj