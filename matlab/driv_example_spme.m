function driv_example_spme
% Example driver for BD simulation.
% The driver program generates the particle configuration and calls bd_driver.
% This version uses SPME.

path(path,'../mex');
rng('default')

npos = 100;                  % number of particles
phi = 0.3;                   % volume fraction
L = (4/3*pi*npos/phi)^(1/3); % box width
pos = rand(3,npos)*L;        % initial positions in box
radii = ones(1,npos);        % particle radii

tic
bd_driver(pos, radii, L, 'test.xyz');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bd_driver(pos, radii, L, filename)

deltat = 0.002;              % time step size
num_steps = 100;             % output interval
num_intervals = 100;         % number of intervals

% repulsion potential is k0/2 (r-r0)^2
r0 = 2;   % A
k0 = 125; % kcal/mol/A^2

% SPME params
xi = 1.5*sqrt(pi)/L;
mode = 'full';
rmax = 5.;
dim = 32;
porder = 4;

spmefun  = @(pos,f) rpy_spme_mex(pos, radii, L, f, mode, xi, rmax, dim, porder);
forcefun = @(pos) force_steric_mex(pos, radii, L, r0, k0);

for i=1:num_intervals
  pos = bd_spme(pos, L, num_steps, deltat, spmefun, forcefun);
  write_xyz(filename, num2str(i*num_steps*deltat), pos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
