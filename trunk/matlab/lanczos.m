function [v h] = lanczos(a, v1, m)
% [v h] = lanczos(a, v1, m)
%
%  Lanczos method for computing a basis for a Krylov subspace, 
%  modified Gram-Schmidt version
%
%  a  = matrix or a function handle
%  v1 = initial vector (normalized v1 will be first vector in basis)
%  m  = size of basis
%  v  = Lanczos basis
%  h  = m by m upper-Hessenberg matrix

%  Ref: Saad, p. 174

if isa(a, 'function_handle')
  afun = a;
else % a is a matrix
  afun = @(x) a*x;
end

n = length(v1);
v = zeros(n,m);
v(:,1) = v1/norm(v1);
h = zeros(m,m);

for j = 1:m
  w = afun(v(:,j));
  if (j > 1)
    w = w - h(j-1,j)*v(:,j-1);
  end

  h(j,j) = w'*v(:,j);
  w = w - h(j,j)*v(:,j);

  if (j+1 <= m)
    h(j+1,j) = norm(w);
    h(j,j+1) = h(j+1,j);
    if (h(j+1,j) < 1e-14)
      fprintf('arnoldi: at step %d, norm of w is %e\n', j, h(j+1,j));
      return
    end
    v(:,j+1) = w/norm(w);
  end
end

% test: the following should be nearly zero:
% norm(v'*afun(v) - h, 'fro')
