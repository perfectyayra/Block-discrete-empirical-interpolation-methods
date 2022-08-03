function irow = cur_deim_block_maxvol(U,k,p,tol)
if nargin < 3 || isempty(p), p = 2; end
if nargin < 4 || isempty(tol), tol = 2e-2; end


irow=zeros(1,k);
for j = 1:(k/p)
    irow((j-1)*p+1:p*j) = maxvol(U(:,(j-1)*p+1:p*j),tol);
	if p*j+p<=k
      U(:,p*j+1:p*j+p) = U(:,p*j+1:p*j+p) - U(:,1:p*j) * (U(irow(1:p*j),1:p*j) \ U(irow(1:p*j),p*j+1:p*j+p));
	end
end



