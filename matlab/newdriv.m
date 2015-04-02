function newdriv(phi, deltat, chol_update, trunc)
% function newdriv(phi, deltat, chol_update)
% phi = 0.3
% deltat = 0.001 or 0.002
% chol_update = 1, 10, 50, 100

path(path,'../mex');
rng('default')

npos = 100;                  % number of particles
phi = 0.3;                   % volume fraction
L = (4/3*pi*npos/phi)^(1/3); % box width
pos = rand(3,npos)*L;        % initial positions in box
radii = ones(1,npos);        % particle radii

%filename = sprintf('p%1dn%03dt%1dc%03d.xyz', floor(10*phi), npos, floor(deltat*1000), chol_update);
filename = sprintf('newtrunc%02dp%1dn%03dt%1dc%03d.xyz', trunc, floor(10*phi), npos, floor(deltat*1000), chol_update);
filename
tic
bd_driver(pos, radii, L, filename, deltat, chol_update, trunc);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = bd_driver(pos, radii, L, filename, deltat, chol_update, trunc)

% All cases take same number of timesteps, even if timestep is smaller.
% But num outputs will be different since output interval is the same.
% Output interval must be at least 100 steps since max chol_update is 100.
% Below gives output interval of 0.2 ps.

if deltat == 0.001
  num_steps = 200;             % output interval
  num_intervals = 100000;      % number of intervals
elseif deltat == 0.002
  num_steps = 100;             % output interval
  num_intervals = 200000;      % number of intervals
elseif deltat == 0.0001
  num_steps = 100;             % output interval
  num_intervals = 200000;      % number of intervals
else
  error('bad deltat');
end

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
  if mod(i,10000) == 0, fprintf('%d of %d\n', i, num_intervals); end
  %pos = mybd(pos, L, num_steps, deltat, matrixfun, forcefun, chol_update);
  pos = mybd_trunc(pos, L, num_steps, deltat, matrixfun, forcefun, chol_update, trunc);
  write_xyz(filename, num2str(i*num_steps*deltat), pos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = rpy_ewald_mex_wrapper(pos, L, xi, nr, nk)
%[mob mob_self mob_real mob_recip] = rpy_ewald(pos, L, xi, nr, nk);
d = rpy_ewald_mex(pos, L, xi, nr, nk);
d = d + d' + rpy_overlap_correction_mex(pos, L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = mybd(pos, L, num_steps, deltat, matrixfun, forcefun, chol_update)
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

cholesky_update_interval = chol_update;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = mybd_trunc(pos, L, num_steps, deltat, matrixfun, forcefun, chol_update, trunc)
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

cholesky_update_interval = chol_update;
% num_steps should be multiple of cholesky update interval

num_outer = num_steps/cholesky_update_interval;

% loop over time steps
for step = 1:num_outer

  pos0 = mod(pos,L);
  d = matrixfun(pos0);

  omega = randn(3*npos,trunc);
  v = mycgs(d*(d*(d*omega)));
  b = v'*d*v;
  [uu ss ~] = svd(b);
  u = v*uu;
  s = sqrtm(ss);
  % no longer need uu and ss

  for inner = 1:cholesky_update_interval
    % compute forces
    forces = forcefun(pos);

    % randomized approximation with good correction
    % implements equation 8 in the proposal
    % u*s*u' = v*sqrtm(b)*v' but the left side is faster because s is diagonal
    z = randn(3*npos,1);
    first_term = u*(s*(u'*z));
    second_term = z - v*(v'*z);  % use v not u
    % possible that z'*d*z changes slowly
    y = augment_vector(first_term, second_term, z'*d*z);

    % update positions (do not apply minimum image convention)
    pos(:) = pos(:) + deltat*(d*forces(:)) + sqrt(2*deltat)*y;
  end

end

function y = augment_vector(v, w, beta2)
% 
a = w'*w;
b = 2*v'*w;
c = v'*v - beta2;
% roots have same magnitude and opposite signs
alpha = (-b + sqrt(b*b-4*a*c))/(2*a);
y = v + alpha*w;

function V = mycgs(A)
% classical Gram schmidt
% only return columns of V
n = size(A,2); % number of cols
for j=1:n
  vj = A(:,j);
  if (j > 1)
     r = V(:,1:j-1)' * A(:,j);  % matvec
     for i = 1:j-1
        vj = vj - r(i) * V(:,i);
     end
  end
  V(:,j) = vj / norm(vj);
end
