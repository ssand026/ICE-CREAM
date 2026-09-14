function [A] = parts(A,indices)
% PARTS(A,ind1,ind2,...) returns an inline index into the array
%
% Arguments after the first correspond to the index into each dimension 
% of the array, where the (n+1)th arguement extracts elements along the
% nth dimension of A. To extract all parts of A along a given dimension, 
% pass the argument "~" or ":", i.e. A(:,1:3,:,5) = parts(A,":",1:3,"~",5)
%
% Additionally, index arguments can be given negative values, where the 
% (k-1)th element from the end corresponds to an index value of -k,
% i.e. A(:,[end-2:end]) = parts(A,":",[-3:-1]).
%
% SEE ALSO: PART
arguments
	A {mustBeNumericOrLogical}
end
arguments (Repeating)
	indices {mustBeA(indices,["numeric","string"]),mustBeNonempty}
end

numInd = numel(indices);

% check number of indices
if isscalar(A) || isempty(A)
	% indexing not possible
	return
elseif ndims(A) == numInd
	% maintain size of A and extract indices
elseif sum(size(A)~=1) == numInd
	% remove dims of size 1 from A and extract indices
	A = squeeze(A);
else
	error("ERROR: invalid number of indices for the array")
end

% check indices for each dimension
dimList = [1:numInd];
for ii = 1:numInd
	x = indices{ii};
	if (isStringScalar(x) && strcmp(x,":")) || (isnumeric(x) && isequal(x,[1:size(A,ii)]))
		% ignore elements along this dimension
		dimList(dimList==ii) = [];
	elseif isnumeric(x) && all(isreal(x)) && all(isfinite(x)) && all(x~=0)
		% extract elements along this dimension
	elseif isstring(x)
		error("ERROR: invalid input string for the array indices")
	elseif isnumeric(x)
		error("ERROR: invalid numeric input for the array indices")
	end
end

% extract indices along each valid dimension
for ii = dimList
	A = part(A,indices{ii},ii);
end
% DONE
end