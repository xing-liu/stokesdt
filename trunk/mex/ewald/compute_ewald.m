function [mat] = compute_ewald(pos, rdi, L, tol, mode, varargin)
% [mat] = compute_ewald(pos, rdi, L, tol, mode, ...)
%
% INPUT
%   pos     the array of particle coordinates
%   rdi     the array of particle radii
%   L       the simulation box size
%   tol     the requested tolerance of Ewald errors
%   mode    'full', 'real' or 'recip'
%   xi      the Ewald parameter (optional);
%           if not specified, the optimal xi will be
%           automatically chosen
%
% OUTPUT
%   mat     the array of the mobility matrix
