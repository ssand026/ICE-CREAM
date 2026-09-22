function [out] = curl(obj)
% computes the curl of a vector field
arguments
	obj FE_coeff
end

% check if the coefficient is a vector
if obj.rank~=1
	error("ERROR: coefficient must be a vector")
elseif size(obj.vals,2)~=1
	% transpose to column-vector
	obj = obj.';
end

% get number of polynomial terms in the gradient
[numPoly,numVtx] = size(obj.poly);

% apply the chain rule: for each polynomial term, take the derivative
% with respect to a single barycoordinate (whose direction is given by bvec)
bvecIndx = [];
polyIndx = [];
dataIndx = [];

for ii = 1:numPoly
	% degree of the curent polynomial term
	deg = nnz(obj.poly(ii,:));
	
	% get the set of valid derivative positions
	row = (1:deg).';
	col = find(obj.poly(ii,:)).';
	vec = accumarray([row,col],1,[deg,numVtx]);
	
	% remaining polynomial powers
	rem = obj.poly(ii,:)-vec;
	polyIndx = [polyIndx; rem]; %#ok
	bvecIndx = [bvecIndx; col]; %#ok
	dataIndx = [dataIndx; repmat(ii,deg,1)]; %#ok
end

% compute the curl coefficients
df = obj.mesh.baryvec(:,:,bvecIndx,:);
if obj.isConstant
	f = obj.vals;
else
	f = obj.vals(:,:,dataIndx,:);
end
data = cross(f,df,1);

% return as a new coefficient
out = FE_coeff(obj.mesh,data,polyIndx,dataIndx);
% DONE
end