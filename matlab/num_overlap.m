function num = num_overlap(pos, rdi, L)
% num = num_overlap(pos, rdi, L)
%  Number of overlapping particles.
%  Warning: this may not work properly for polydisperse cases.

npos = size(pos,2);
[ia ja] = pairlist_mex(pos, rdi, L, 2, 'R');
% returns overlapping pairs (including a particle with itself)
num = length(ja) - npos;
