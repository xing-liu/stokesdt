function forces = force_repulsion_parallel(pos, L, krepulsion)
% forces = force_repulsion_parallel(pos, L, krepulsion)
%  pos should be 3 by npos array

npos = size(pos,2);
forces = zeros(3,npos);
radii = ones(npos,1);
cutoff = 2;
cutoff2 = cutoff*cutoff;

[ia ja] = pairlist_mex(pos, radii, L, cutoff, 'R');
% 0-based indexing

length(ja)

for i=1:npos
  for j=ia(i)+1:ia(i+1)
    r = pos(:,i) - pos(:,ja(j)+1);

    % compute minimum image
    r = mod(r, L);  % assumes mod(-0.1,1)=0.9
    if r(1) > L/2, r(1) = r(1)-L; end
    if r(2) > L/2, r(2) = r(2)-L; end
    if r(3) > L/2, r(3) = r(3)-L; end

    s2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3);

    if (s2 == 0)
      continue;
    end

    if (s2 <= cutoff2)
      s = sqrt(s2);
      e = r/s;  % unit vector from j to i
      f = krepulsion*(2 - s);
      forces(:,i) = forces(:,i) + f*e;
    end

  end
end
