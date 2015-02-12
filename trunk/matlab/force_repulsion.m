function forces = force_repulsion(pos, L, krepulsion)
% forces = force_repulsion(pos, L, krepulsion)

npos = size(pos,2);
forces = zeros(3,npos);

for i = 1:npos
  posi = pos(:,i);
  for j = i+1:npos
    posj = pos(:,j);
    r = posi-posj;

    % compute minimum image
    r = mod(r, L);  % assumes mod(-0.1,1)=0.9
    if r(1) > L/2, r(1) = r(1)-L; end
    if r(2) > L/2, r(2) = r(2)-L; end
    if r(3) > L/2, r(3) = r(3)-L; end

    s2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3);

    if (s2 < 4)
      s = sqrt(s2);
      e = r/s;  % unit vector from j to i
      f = krepulsion*(2 - s);
      forces(:,i) = forces(:,i) + f*e;
      forces(:,j) = forces(:,j) - f*e;
    end
  end
end
