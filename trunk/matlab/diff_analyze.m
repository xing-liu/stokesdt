function diff_analyze(dvec, numseg)
n = length(dvec);
seglen = floor(n/numseg);
s = zeros(numseg,1);
for i=1:numseg
  s(i) = mean(dvec((i-1)*seglen+1:i*seglen));
end
fprintf('Original number of samples: %d\n', n);
fprintf('Discarding first %d samples\n', seglen);
averg = mean(s(2:end));
stdev = std(s(2:end));
fprintf('Mean: %f    std: %f\n', averg, stdev);
fprintf('Error bars:  %f  %f\n', averg-stdev, averg+stdev);
s
