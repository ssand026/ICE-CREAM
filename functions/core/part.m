function [A] = part(A,indx,dim)
% PART(A,ind,dim) extracts the specified parts of A along the given dimension
% SEE ALSO: PARTS
arguments
	A {mustBeNumericOrLogical,mustBeNonempty}
	indx (:,1) {mustBeInteger,mustBeNonzero,mustBeFinite,mustBeNonempty}
	dim (1,1)  {mustBeInteger,mustBePositive}
end

if dim > ndims(A)
	error("ERROR: specified dim is larger than the number of array dimensions")
elseif max(indx) > +size(A,dim)
	error("ERROR: positive value/s in the index exceed the dimension length")
elseif min(indx) < -size(A,dim)
	error("ERROR: negative value/s in the index exceed the dimension length")
end
% convert negative indices into distance from the end, i.e. mat(:,-2) = mat(:,end-1); 
indx(indx < 0) = indx(indx < 0) + size(A,dim) + 1;

% permute array so the indexed-over dimension is the first one, extract the
% relevant indices, reshape and un-permute
fwd = circshift([1:ndims(A)],-(dim-1));
bwd = circshift([1:ndims(A)],+(dim-1));

matSize = size(A);
A = permute(A,fwd);
A = reshape(A(indx,:),[numel(indx),matSize(fwd(2:end))]);
A = permute(A,bwd);
end