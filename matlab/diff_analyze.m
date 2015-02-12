function diff_analyze(dvec)
% assuming 9999 data points
% throw away first part and divide remaining into three segments
s(1) = mean(dvec(4000:5999));
s(2) = mean(dvec(6000:7999));
s(3) = mean(dvec(8000:9999));
fprintf('mean: %f  std: %f\n', mean(s), std(s));
