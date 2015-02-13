function pos = bd_spme(pos, L, num_steps, deltat, spmefun, forcefun)
% pos = bd(pos, L, num_steps, deltat, matrixfun, forcefun)
%   Brownian dynamics simulation with periodic boundary conditions.
%
% pos       = particle positions (3 x npos) (Angstroms)
% L         = cubical box side length (Angstroms)
% num_steps = number of simulation steps to perform
% deltat    = time step size (ps)
% spmefun   = function returning velocities given forces
% forcefun  = function returning external forces, given positions

% number of particles
npos = size(pos,2);
if size(pos,1) ~= 3
  error('leading dimension of pos must be 3')
end

% loop over time steps
for step = 1:num_steps
  step

  % not known if needed
  pos0 = mod(pos,L);

  % Brownian displacement
  % does spme need pos def correction?
  z = randn(3*npos,1);
  [v h] = lanczos(@(f)spmefun(pos0,f), z, 5); % hardcode size of basis
  sqrth = sqrtm(h);
  if ~isreal(sqrth(:,1))
    eig(h)
    error('sqrth is not real'); 
  end
  y = v*sqrth(:,1)*norm(z);

  % compute velocity due to external forces
  forces = forcefun(pos);
  dtimesf = spmefun(pos, forces(:));
  if ~isreal(dtimesf), error('dtimesf is not real'); end

  % update positions (do not apply minimum image convention)
  pos(:) = pos(:) + deltat*dtimesf + sqrt(2*deltat)*y;

end

%%%%%

