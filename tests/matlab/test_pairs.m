clear all
path(path,'../../install/mex');
rng('default')

% L x L x L simulation box
L = 1;

% number of particles
npos = 100;

% random initial positions in simulation box
pos = rand(3,npos)*L;
radii = ones(npos,1);

% interaction cutoff
cutoff = 0.3;
% cutoff cannot be greater than L/2

[ia ja] = compute_pairlist(pos(:), radii, L, cutoff, 'R');

% change to 1-based indexing
ia = ia + 1;
ja = ja + 1;

mark = zeros(npos,npos);
for i=1:npos
  for j=ia(i):ia(i+1)-1
    mark(i,ja(j)) = 1;
  end
end

% note: both i,j and j,i are in pairlist
% note: diagonal entry is stored
if (nnz(mark) ~= length(ja))
  error('duplicate entry in pairlist');
end

% construct pairlist using brute force
% must compute distances using nearest image
brute = zeros(npos,npos);
for j=1:npos
  for i=j:npos
    r = pos(:,i) - pos(:,j);

    % compute minimum image
    r = mod(r, L);  % assumes mod(-0.1,1)=0.9
    if r(1) > L/2, r(1) = r(1)-L; end
    if r(2) > L/2, r(2) = r(2)-L; end
    if r(3) > L/2, r(3) = r(3)-L; end

    s2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3);

    if (s2 <= cutoff*cutoff)
      brute(i,j) = 1;
    end
  end
end
brute = spones(brute+brute');

% compare
if (norm(mark-brute,'fro') == 0)
  fprintf('PASS\n');
else
  error('FAIL');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test polydisperse case
