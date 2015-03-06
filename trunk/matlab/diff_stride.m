function diffstride = diff_stride(filename, maxstride, first)
% diffstride = diff_stride(filename, maxstride, first)
% compute diffusion coefficient for different time intervals (strides)
% last two arguments are optional

if nargin < 3
  % first frame to analyze (skip earlier ones)
  first = 30000;
end
if nargin < 2
  maxstride = 3000;
end

% read number of particles and time between frames
fid = fopen(filename, 'r');
[npos        count] = fscanf(fid, '%d', 1);
[framedeltat count] = fscanf(fid, '%f', 1);
fclose(fid);

% read number of frames
[status result] = system(sprintf('wc %s', filename));
numlines = sscanf(result, '%d', 1);
numframes = floor(numlines/(npos+2));

fprintf('num particles: %d\n', npos);
fprintf('framedeltat: %f\n', framedeltat);
fprintf('num frames: %d\n', numframes);
fprintf('firstframe: %d\n', first);

all = zeros(npos,3,numframes);
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
  diffstride(stride) =  tot/(6*framedeltat*stride*num);
  if (mod(stride,100)==0)
    fprintf('stride %4d   %f\n', stride, diffstride(stride));
  end
end
