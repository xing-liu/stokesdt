function diffstride = diff_stride
% compute diffusion coefficient for different time intervals (strides)

maxstride = 300;
numframes = 42400;
filename = 'phi30-0100.chol1.xyz';
filename = 'notrunc-0100.xyz';
filename = 'phi30-0100.halftimestep.xyz';
deltat = 0.1; % from input file
first = 10000; % first frame to analyze (skip earlier ones)

% set parameters above this line

all = zeros(100,3,numframes);
fid = fopen(filename, 'r');
for i=1:numframes
  if mod(i,1000)==0, fprintf('%d ', i); end
  [npos count] = fscanf(fid, '%d', 1);        if (count < 1),      break, end
  [time count] = fscanf(fid, '%f', 1);        if (count < 1),      break, end
  [pos  count] = fscanf(fid, '%f', [4 npos]); if (count < 4*npos), break, end
  all(:,:,i) = pos(2:4,:)';
end
fclose(fid);
fprintf('\n');

diffstride = zeros(maxstride,1);
for stride=1:maxstride
  pos0 = all(:,:,first);
  tot = 0;
  num = 0;
  for i=first+stride:stride:numframes
    pos = all(:,:,i);
    r = pos - pos0;
    s = r(:,1).*r(:,1) + r(:,2).*r(:,2) + r(:,3).*r(:,3);
    tot = tot + mean(s);
    num = num + 1;
    pos0 = pos;
  end
  diffstride(stride) =  tot/(6*deltat*stride*num);
  fprintf('stride %3d   %f\n', stride, diffstride(stride));
end
