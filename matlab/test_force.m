clear all
path(path,'../mex');
rng('default')

npos = 100;
phi = .3;
particle_vol = npos*4/3*pi;
L = (1/phi)*particle_vol^(1/3);
krepulsion = 125;

pos = rand(3,npos)*L;

tic; f1 = force_steric_mex(pos, ones(npos,1), L, 2., 125.); toc
tic; f2 = force_repulsion(pos, L, 125.); toc
tic; f3 = force_repulsion_mex(pos, L, 125.); toc
norm(f1(:)-f2(:))
norm(f3(:)-f2(:))

% up to 200 particles, force_repulsion_mex is fastest
% omp version of this has high overhead - may not pay off for larger problems
% before the verlet list version is faster
