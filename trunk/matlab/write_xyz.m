function write_xyz(filename, label, pos)
% write_xyz(filename, label, pos)

npos = size(pos, 2);
fid = fopen(filename, 'a'); % file opened each time
fprintf(fid, '%d\n', npos);
fprintf(fid, '%s\n', label);
for i=1:npos
  fprintf(fid, '%s %f %f %f\n', 'X', pos(:,i));
end
fclose(fid);
