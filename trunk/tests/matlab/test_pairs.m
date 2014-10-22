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
cutoff = 0.1;

% pi is not a good name
[indi indj] = compute_verlet(pos(:), radii, L, cutoff, 'R');

% change to 1-based indexing
indi = indi + 1;
indj = indj + 1;

% mark the pairs in a lower triangular array
mark = zeros(npos,npos);
for i=1:length(indi)
  mark(indi(i),indj(i)) = 1;
end

% how are the pairs ordered?  Lower triang?  Upper triang?
% can they appear more than once?  will i,j and j,i both be output?
%%% all pairs output...

% check answer by computing distance between all pairs
% must compute distances using nearest image
for j=1:npos
  for i=j+1:npos
    d = norm(pos(:,indi(i)) - pos(:,indj(i)));
    m = mark(indi(i),indj(i));
    if (d <= cutoff && m ~= 1)
      error('test failure 1');
    end
    if (d > cutoff && m ~= 0)
      d
      m
      error('test failure 2');
    end
  end
end

% test polydisperse case
