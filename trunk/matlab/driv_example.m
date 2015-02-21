function driv_example(phi)
% Example driver for BD simulation.
% The driver program generates the particle configuration and calls bd_driver.

path(path,'../mex');
rng('default')

npos = 100;                  % number of particles
phi = 0.3;                   % volume fraction
L = (4/3*pi*npos/phi)^(1/3); % box width
pos = rand(3,npos)*L;        % initial positions in box
radii = ones(1,npos);        % particle radii

filename = sprintf('phi%2d-0100.xyz', floor(100*phi));
tic
bd_driver(pos, radii, L, 'phi30-0100.xyz');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = bd_driver(pos, radii, L, filename)

deltat = 0.002;              % time step size
num_steps = 1000;            % output interval
num_intervals = 1000000;     % number of intervals

% Ewald params
xi = 1.5*sqrt(pi)/L;
nr = 2;
nk = 3;

% repulsion potential is k0/2 (r-r0)^2
r0 = 2;   % A
k0 = 125; % kcal/mol/A^2

matrixfun = @(pos) rpy_ewald_matrix_mex(pos, radii, L, 1e-4, 'full', xi);
matrixfun = @(pos) rpy_ewald_mex_wrapper(pos, L, xi, nr, nk);
forcefun  = @(pos) force_repulsion(pos, L, k0); % very slow
forcefun  = @(pos) force_steric_mex(pos, radii, L, r0, k0);
forcefun  = @(pos) force_repulsion_mex(pos, L, k0);

for i=1:num_intervals
  pos = mybd(pos, L, num_steps, deltat, matrixfun, forcefun);
  write_xyz(filename, num2str(i*num_steps*deltat), pos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = rpy_ewald_mex_wrapper(pos, L, xi, nr, nk)
%[mob mob_self mob_real mob_recip] = rpy_ewald(pos, L, xi, nr, nk);
d = rpy_ewald_mex(pos, L, xi, nr, nk);
d = d + d' + rpy_overlap_correction_mex(pos, L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = mybd(pos, L, num_steps, deltat, matrixfun, forcefun)
% pos = bd(pos, L, num_steps, deltat, matrixfun, forcefun)
%   Brownian dynamics simulation with periodic boundary conditions.
%
% pos       = particle positions (3 x npos) (Angstroms)
% L         = cubical box side length (Angstroms)
% num_steps = number of simulation steps to perform
% deltat    = time step size (ps)
% matrixfun = function returning hydrodynamic mobility matrix, given positions
% forcefun  = function returning external forces, given positions

% number of particles
npos = size(pos,2);
if size(pos,1) ~= 3
  error('leading dimension of pos must be 3')
end

cholesky_update_interval = 50;
% num_steps should be multiple of cholesky update interval

num_outer = num_steps/cholesky_update_interval;

% loop over time steps
for step = 1:num_outer

  pos0 = mod(pos,L);
  d = matrixfun(pos0);

  chold = chol(d)';

  for inner = 1:cholesky_update_interval
    % compute forces
    forces = forcefun(pos);

    % Brownian displacement
    z = randn(3*npos,1);
    y = chold*z;

    % update positions (do not apply minimum image convention)
    pos(:) = pos(:) + deltat*(d*forces(:)) + sqrt(2*deltat)*y;
  end

end
