path(path,'../../install/mex');
path(path,'/home/edmond/bd-2014/mex');

npos = 1000;
phi = .3;
particle_vol = npos*4/3*pi;
L = (1/phi)*particle_vol^(1/3);
krepulsion = 125;

pos = rand(3,npos)*L;

forces1 = force_repulsion(pos, L, krepulsion);
forces2 = force_repulsion_parallel(pos, L, krepulsion);

norm(forces1-forces2)

