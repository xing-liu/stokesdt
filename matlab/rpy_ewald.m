function [mob mob_self mob_real mob_recip] = rpy_ewald(pos, L, xi, nr, nk)
% mob = rpy_ewald(pos, L, xi, nr, nk)
%   RPY mobility matrix with Ewald sum.

% dimensionless particle radii
a = 1;
a3 = 1;

% wrap positions into box
pos = mod(pos, L);
npos = size(pos,2);

% dense mobility matrix
mob_real  = zeros(npos*3);
mob_recip = zeros(npos*3);
eye3 = eye(3);

% self-term
mob_self = (1 - (6*a - 40/3*xi*xi*a3)*xi/sqrt(pi)) * eye(npos*3);

% real-space sum
% loop over particle pairs including with self
for i = 1:npos
  posi = pos(:,i);
  for j = i:npos
    posj = pos(:,j);
    rvec0 = posi - posj;

    temp = zeros(3);

    for x = -nr:nr
    for y = -nr:nr
    for z = -nr:nr

      rvec = rvec0 + [x y z]'*L;

      r = norm(rvec);
      if (r == 0)
        continue; % own particle
      end
      e = rvec/r;
          
      % sum tensors into matrix
      [m11 m12] = scalar_rpy_ewald_real(r, xi);
      temp = temp + m11*eye3 + m12*e*e';

    end
    end
    end

    if (i == j)
      temp = temp/2;
    end

    id = i*3-2:i*3;
    jd = j*3-2:j*3;

    mob_real(id,jd) = temp;

  end
end

% reciprocal-space sum
% loop over particle pairs including self
for i = 1:npos
  posi = pos(:,i);
  for j = i:npos

    posj = pos(:,j);
    rvec = posi - posj;

    temp = zeros(3);

    for x = -nk:nk
    for y = -nk:nk
    for z = -nk:nk

      kvec = 2*pi/L*[x y z]';
      k = norm(kvec);
      if (k == 0)
        continue;
      end
      e = kvec/k;
          
      % sum tensors into matrix
      % scalars at each recip lattice point can be saved
      m2 = scalar_rpy_ewald_recip(k, xi);
      temp = temp + m2*cos(kvec'*rvec)*(eye3 - e*e');

    end
    end
    end

    id = i*3-2:i*3;
    jd = j*3-2:j*3;

    if (i == j)
      temp = temp/2;
    end

    mob_recip(id,jd) = temp;

  end
end

mob_recip = mob_recip / (L*L*L);
mob_recip = mob_recip + mob_recip';
mob_real  = mob_real  + mob_real';

mob = mob_real + mob_recip + mob_self;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m11 m12] = scalar_rpy_ewald_real(r, xi) 
% [m11 m12] = scalar_rpy_ewald_real(r, xi)

a = 1;
a3 = 1;

xi2 = xi*xi;
xi3 = xi2*xi;
xi5 = xi3*xi2;
xi7 = xi5*xi2;

r2 = r*r;
r4 = r2*r2;
ri = 1/r;
ri2 = ri*ri;
ri3 = ri*ri2;

erfc_xi_r = erfc(xi*r);
pi_exp = 1/sqrt(pi) * exp(-xi2*r2);

m11 = (0.75*a*ri + 0.5*a3*ri3)*erfc_xi_r + ( 4*xi7*a3*r4 + 3*xi3*a*r2 - 20*xi5*a3*r2 - 4.5*xi*a + 14*xi3*a3 +   xi*a3*ri2)*pi_exp;
m12 = (0.75*a*ri - 1.5*a3*ri3)*erfc_xi_r + (-4*xi7*a3*r4 - 3*xi3*a*r2 + 16*xi5*a3*r2 + 1.5*xi*a -  2*xi3*a3 - 3*xi*a3*ri2)*pi_exp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m2 = scalar_rpy_ewald_recip(k, xi)
% m2 = scalar_rpy_ewald_recip(k, xi)

a = 1;
a3 = 1;

k2 = k*k;
xii2k2 = k2/(xi*xi);

m2 = (a-a3*k2/3) * (1 + 0.25*xii2k2 + 0.125*xii2k2*xii2k2) * 6*pi/k2 * exp(-0.25*xii2k2);
