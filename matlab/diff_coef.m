function dvec = diff_coef(filename)
% dvec = diff_coef(filename)
%  Compute vector of diffusion coefficients (one for each frame except first)
%  from xyz file.  Input xyz file must have physical time as the label.
%  This function is designed to work on a partially completed xyz file
%  that is still being used as output from a simulation.

% If label does not contain physical time, then code could be modified
% to use the physical time between frames.

% preallocate space for dvec
dvec = zeros(1000000,1);

fid = fopen(filename, 'r');
framenum = 0;
while true
  [npos count] = fscanf(fid, '%d', 1);        if (count < 1),      break, end
  [time count] = fscanf(fid, '%f', 1);        if (count < 1),      break, end
  [pos  count] = fscanf(fid, '%f', [4 npos]); if (count < 4*npos), break, end
  framenum = framenum + 1;
  if (framenum == 1)
    pos0 = pos;
    time0 = time;
  else
    r = (pos - pos0)';
    % first column of r is particle symbol
    s = r(:,2).*r(:,2) + r(:,3).*r(:,3) + r(:,4).*r(:,4);
    dvec(framenum-1) = mean(s)/(6*(time-time0));
    pos0 = pos;
    time0 = time;
  end
end
fclose(fid);
if (framenum < 1)
  error('Need more than one frame in xyz file');
end
dvec = dvec(1:framenum-1);
fprintf('Number of frames: %d\n', framenum);
