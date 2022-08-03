function irow = block_deim_qr(U, k, p)
% U is singular vectors
% k desired number of indices
% p block size
% Revision date: June 18, 2022
% (C) Perfect Gidisu, Michiel Hochstenbach 2022

if nargin < 3|| isempty(p), p = 2; end


irow=zeros(1,k);
for j = 1:round(k/p)
    irow((j-1)*p+1:p*j)=qdeim(U(:,(j-1)*p+1:p*j),p);
	if p*j+p<=k
      U(:,p*j+1:p*j+p) = U(:,p*j+1:p*j+p) - U(:,1:p*j) * (U(irow(1:p*j),1:p*j) \ U(irow(1:p*j),p*j+1:p*j+p));
	end 
    
end
