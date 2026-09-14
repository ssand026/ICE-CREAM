function [vals,poly,index] = collect(vals,poly,index)
% combine coefficient terms with the same polynomial components
arguments
	vals (:,:,:,:) double
	poly (:,:) double
	index (:,:) double = [];
end

% circumvent the costly overhead of "unique(...)" using:
%	[C,ia,ic] = matlab.internal.math.uniquehelper(A,doSort,isFirst,byRows)
fast_unique = @(x,opts) matlab.internal.math.uniquehelper(x,opts(1),opts(2),opts(3));

% find identical polynomial terms
index_poly = full([index,poly]);
[index_poly,~,ordering] = fast_unique(index_poly,[false,true,true]);

% check if there are any like terms to combine
originalLen = numel(ordering);
combinedLen = max(ordering);

if (combinedLen < originalLen)
	% split the reduced index_poly back into the individual components
	numIndex = size(index,2);
	index = index_poly(:,1:numIndex);
	poly = index_poly(:,numIndex+1:end);
	
	% push the polynomial dim to the front and flatten
	dimOrder = [3,1,2,4:ndims(vals)];
	valsSize = size(vals,dimOrder(2:end));
	
	vals = permute(vals,dimOrder);
	vals = reshape(vals,originalLen,[]);

	% generate a matrix to combine the identical polynomial terms
	combMat = sparse(ordering,(1:originalLen),true,combinedLen,originalLen);
	vals = combMat * vals;

	% unflatten the value array and undo the permutation
	vals = reshape(vals,[combinedLen,valsSize]);
	[~,invOrder] = sort(dimOrder);
	vals = permute(vals,invOrder);
end
% DONE
end

