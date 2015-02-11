function pos = bd(pos, L, num_steps, deltat, matrixfun, forcefun)
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
