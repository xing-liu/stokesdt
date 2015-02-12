function a = rpy_overlap_correction(pos, L)
% a = rpy_overlap_correction(pos, L)
%  a = correction matrix

npos = size(pos,2);
a = zeros(3*npos);
eye3 = eye(3);

for i=1:npos
  posi = pos(:,i);
  for j=i+1:npos
    posj = pos(:,j);
    r = posi-posj;

    % compute minimum image
    r = mod(r, L);  % assumes mod(-0.1,1)=0.9
    if r(1) > L/2, r(1) = r(1)-L; end
    if r(2) > L/2, r(2) = r(2)-L; end
    if r(3) > L/2, r(3) = r(3)-L; end

    s2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3);

    % sum tensors into matrix
    if (s2 < 4)
      s = sqrt(s2);
      e = r/s;

      id = i*3-2:i*3;
      jd = j*3-2:j*3;

      t = 0.09375*s; % 3/32*s

      a(id,jd) = ( (1-3*t) - 0.75/s*(1+2/(3*s2)) ) * eye3 ...
               + (      t  - 0.75/s*(1-2/s2)     ) * e*e';
    end
  end
end

a = a + a';
