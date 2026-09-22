function [out] = isunique(x,dim)
% ISUNIQUE Checks if slices of an array along the given dimension only occur once
% 
% SYNTAX:
%	isunique(x) OR isunique(x,dim) where (dim < 1)
%		checks if every element of 'x' is unique, output is the same size as 'x'
%	isunique(x,dim) where (dim >= 1)
%		slices 'x' along the specified dimension,
%	isunique(A,[])
%		checks slices of 'x' along the first dimension with non-unit length.
%
arguments
	x {mustBeA(x,{'numeric','logical','string'})}
	dim uint32 {mustBeScalarOrEmpty,mustBeInteger} = 0;
end

if isempty(x)
	% no unique elements
	out = [];
	return
elseif isscalar(x)
	% unique by definition
	out = true;
	return
elseif isempty(dim)
	% set dim to the first non-singular dimension
	dim = find(size(x)>1,1,"first");
end

if (dim == 0)
	% check uniqueness of every element
	out = isunique(x(:),1);
	out = reshape(out,size(x));
	return
elseif size(x,dim) == 1
	% unique by definition
	out = true;
	return

end

% restructure the array so that each row contains a list of the
% elements found in each slice
numSlices = size(x,dim);
x = permute(x,[1:(dim-1),(dim+1):ndims(x),dim]);
x = permute(reshape(x,[],numSlices),[2,1]);

% sort rows and check if adjacent rows are equal
[sortedA,order] = sortrows(x);
out = any(sortedA(1:end-1,:) ~= sortedA(2:end,:),2);
out = ([true; out(:)] & [out(:); true]);

% undo the row sorting
[~,invOrder] = sort(order);
out = out(invOrder(:));
% DONE
end