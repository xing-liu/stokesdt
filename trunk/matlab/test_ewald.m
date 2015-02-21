clear all
path(path,'../mex');
rng('default')

npos = 100;
phi = 0.3;
L = (4/3*pi*npos/phi)^(1/3);
pos = rand(3,npos)*L;

xi = 1.5*sqrt(pi)/L;
nr = 2;
nk = 3;

tic
corr = rpy_overlap_correction(pos, L);
toc

tic
corr2 = rpy_overlap_correction_mex(pos, L);
toc

norm(corr-corr2,'fro')

% how are nr and nk chosen?
% how can we find out what was chosen
tic
mob0 = rpy_ewald_matrix_mex(pos, ones(npos,1), L, 1e-4, 'full', xi);
toc
% already has correction
% mode is full, real, or recip
% where is the self interaction added?
% what is the meaning of the tolerance?
% other functions below can only handle monodisperse
% how does it handle polydisperse?
% how is it parallelized?

tic
[mob1 mob_self mob_real mob_recip] = rpy_ewald(pos, L, xi, nr, nk);
toc
mob1 = mob1 + corr;

tic
[mob2 mob_self mob_real mob_recip] = rpy_ewald_refactored2(pos, L, xi, nr, nk);
toc
mob2 = mob2 + corr;

tic
mob3 = rpy_ewald_mex(pos, L, xi, nr, nk);
toc
mob3 = mob3 + mob3';
mob3 = mob3 + corr;

norm(mob1-mob0,'fro')
norm(mob1-mob2,'fro')
norm(mob1-mob3,'fro')
