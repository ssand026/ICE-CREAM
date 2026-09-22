function [obj] = by_element(mesh,elements)
% constructs a FE_coeff from a list of values/expressions corresponding to each element of the 
% coefficient vector/tensor field
arguments
	mesh FE_mesh
	elements (:,:) cell
end

[numPts,numTri,numDim] = mesh.num("pts","tri","dim");
coeffDims = size(elements,[1,2]);

% check the list of elements
if ~ismember(coeffDims,[1,numDim])
	error("ERROR: the list of elements has an incorrect size for a valid tensor")
end

% evaluate tensor components
isDisjoint = false(coeffDims);
isConstant = false(coeffDims);
for ii = 1:coeffDims(1)
	for jj = 1:coeffDims(2)
		data_ij = mesh.evalfun(elements{ii,jj});
		[elements{ii,jj},isDisjoint(ii,jj),isConstant(ii,jj),rank] = check(mesh,data_ij);
		if rank ~= 0
			error("ERROR: element ("+ii+","+jj+") must evaluate to a rank-0 tensor")
		end
	end
end

% choose smallest compatible coefficient form to utilize
if any(isDisjoint(:) & ~isConstant(:))
	% must expand all
	for nn = 1:numel(elements)
		elements{nn} = expand(mesh,elements{nn});
	end
elseif all(isConstant(:))
	% must expand non-disjointed
	for nn = find(~isDisjoint)
		elements{nn} = repmat(elements{nn},1,1,1,numTri);
	end
elseif all(~isDisjoint(:))
	% must expand constants
	for nn = find(isConstant)
		elements{nn} = repmat(elements{nn},1,1,numPts,1);
	end
end

% combine into single coefficient
data = zeros([coeffDims, size(elements{1,1},[3,4])]);
for ii = 1:coeffDims(1)
	for jj = 1:coeffDims(2)
		data(ii,jj,:,:) = elements{ii,jj};
	end
end

% invoke constructor
[obj] = FE_coeff(mesh,data);

% DONE
end