function driv_example
% Example driver for BD simulation.
% The driver program generates the particle configuration and calls bd_driver.

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

% Ewald params
xi = 1.5*sqrt(pi)/L;
nr = 2;
nk = 3;

% repulsion potential is k0/2 (r-r0)^2
r0 = 2;   % A
k0 = 125; % kcal/mol/A^2

matrixfun = @(pos) rpy_ewald_mex_wrapper(pos, L, xi, nr, nk);
matrixfun = @(pos) rpy_ewald_matrix_mex(pos, radii, L, 1e-4, 'full', xi);
forcefun  = @(pos) force_repulsion(pos, L, k0); % very slow
forcefun  = @(pos) force_repulsion_mex(pos, L, k0);
forcefun  = @(pos) force_steric_mex(pos, radii, L, r0, k0);

for i=1:num_intervals
  pos = bd(pos, L, num_steps, deltat, matrixfun, forcefun);
  write_xyz(filename, num2str(i*num_steps*deltat), pos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = rpy_ewald_mex_wrapper(pos, L, xi, nr, nk)
%[mob mob_self mob_real mob_recip] = rpy_ewald(pos, L, xi, nr, nk);
d = rpy_ewald_mex(pos, L, xi, nr, nk);
d = d + d' + rpy_overlap_correction(pos, L);

