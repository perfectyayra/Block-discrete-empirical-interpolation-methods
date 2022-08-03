function irow = cur_adap_blk_qr(U,k,p,rho)
if nargin < 3|| isempty(p), p = 2; end
if nargin < 4|| isempty(rho), rho = 0.9; end


irow = deim_adapt(U, k, p, rho);


function index = deim_adapt(U, k, p, rho)
j = 1; index = zeros(1,k);
while j <= k
  if j > 1, U(:,j) = U(:,j) - U(:,1:j-1)*(U(index(1:j-1),1:j-1) \ U(index(1:j-1),j)); end
  [u, ind] = sort(abs(U(:,j)), 'descend');
  if (j+p-1 > k) || (u(2) < rho*u(1))  % Standard (vector) variant
    index(j) = ind(1); j = j+1;
  else                              % Block variant
    if j > 1, U(:,j+1:j+p-1) = U(:,j+1:j+p-1) - U(:,1:j-1)*(U(index(1:j-1),1:j-1) \ U(index(1:j-1),j+1:j+p-1)); end
    [~,~,P]=qr(U(:,j:j+p-1)','vector');
    index(j:j+p-1) = P(1:p); j = j+p;
  end
end



