function index = adap_blk_maxvol(U,k,p,rho,tol)
% U is singular vectors
% k desired number of indices
% p block size
% rho evaluation criteria for switch between standard and block deim
% tol stopping cirteria for maxvol algorithm 
% see also deim and block_deim_maxvol
% Revision date: June 18, 2022
% (C) Perfect Gidisu, Michiel Hochstenbach 2022

if nargin < 3|| isempty(p), p = 2; end
if nargin < 4|| isempty(rho), rho = 0.95; end
if nargin < 5 || isempty(tol), tol = 1e-2; end

j = 1; index = zeros(1,k);
while j <= k
  if j > 1, U(:,j) = U(:,j) - U(:,1:j-1)*(U(index(1:j-1),1:j-1) \ U(index(1:j-1),j)); end
  [u, ind] = sort(abs(U(:,j)), 'descend');
  if (j+p-1 > k) || (u(2) < rho*u(1))  % Standard (vector) variant
    index(j) = ind(1); j = j+1;
  else                              % Block variant
    if j > 1, U(:,j+1:j+p-1) = U(:,j+1:j+p-1) - U(:,1:j-1)*(U(index(1:j-1),1:j-1) \ U(index(1:j-1),j+1:j+p-1)); end
    index(j:j+p-1) = maxvol(U(:,j:j+p-1),tol); j = j+p;
  end
end



