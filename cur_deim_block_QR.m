function irow = cur_deim_block_QR(U, k, p)

if nargin < 3|| isempty(p), p = 2; end


irow=zeros(1,k);
for j = 1:round(k/p)
    irow((j-1)*p+1:p*j)=qdeim(U(:,(j-1)*p+1:p*j),p);
	if p*j+p<=k
      U(:,p*j+1:p*j+p) = U(:,p*j+1:p*j+p) - U(:,1:p*j) * (U(irow(1:p*j),1:p*j) \ U(irow(1:p*j),p*j+1:p*j+p));
	end 
    
end
