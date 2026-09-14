function [data,isDisjoint,isConstant,rank] = check(mesh,data)
% check - Checks if the given array is a valid mesh coefficient
%
% INPUTS:
%	mesh: FE_mesh object
%	data: numeric array representing a scalar/vector/tensor field
%	strictDims(optional,bool): require coefficient to be a pure field (no extra dims)
%
% OUTPUTS:
%	1) data: the re-ordered coefficient array
%	2) isDisjoint: whether the coefficient is discontinuous
%   3) isConstant: whether the coefficient is contant/variable
%	4) rank: the tensor-order of the coefficient
%		- 0: coeff is a scalar field, [d1,d2]==[1,1]
%		- 1: coeff is a vector field, [d1,d2]==[numDim,1]
%		- 2: coeff is a tensor field, [d1,d2]==[numDim,numDim]
%
% SEE ALSO: FE_COEFF
arguments
	mesh FE_mesh
	data double
end

% get the number of spatial dimensions, nodes, and simplices in the mesh
[numPts,numTri,numDim] = mesh.num("pts","tri","dim");

% check if data is already in the proper order
if (size(data,1)==1 || size(data,1)==numDim) && (size(data,2)==1 || size(data,2)==numDim)
	rank = sum(size(data,[1,2])~=1);
	if ndims(data)==3 && size(data,3)==numPts
		% field is variable over the domain
		isDisjoint = false;
		isConstant = false;
		return
	elseif ndims(data)==4 && size(data,4)==numTri
		isDisjoint = true;
		if size(data,3)==1
			% field is constant over each simplex
			isConstant = true;
		else
			% field is variable over each simplex
			isConstant = false;
		end
		return
	elseif ismatrix(data)
		% field is constant over the domain
		isDisjoint = false;
		isConstant = true;
		if rank==1; data = data(:); end
		return
	end

end

% otherwise attempt to permute the array dimensions
%-----------------------------------------------------------
% ensure that the mesh measures are distinct
if (numDim == numPts) || (numDim == numTri) || (numPts == numTri)
	error("ERROR: overlap between the number of nodes, dimensions, or" + ...
		" simplexes in the mesh. These numbers are typically distinct" + ...
		" outside of extremely small meshes")
elseif (numDim == 1)
	error("ERROR: one-dimensional coefficients are not supported");
end

% get size of data
data = squeeze(data);
where_nDim = find(size(data)==numDim);
where_nPts = find(size(data)==numPts);
where_nTri = find(size(data)==numTri);
where_othr = find(~ismember(size(data),[1,numDim,numPts,numTri]));

% check if data dimensions correspond to a valid field type
shape = [numel(where_nPts), numel(where_nTri), numel(where_othr)];
isDisjoint = true;
isConstant = false;
if isequal(shape,[0,0,0])
	% field is constant over the domain
	isDisjoint = false;
	isConstant = true;
elseif isequal(shape,[1,0,0])
	% field is variable over the domain
	isDisjoint = false;
	isConstant = false;
elseif isequal(shape,[0,1,0])
	% field is constant over each simplex
	isDisjoint = true;
	isConstant = true;
elseif isequal(shape,[0,1,1])
	% field is variable over each simplex
	isDisjoint = true;
	isConstant = false;
elseif isequal(shape,[1,0,1]) || (shape(3) > 1)
	% coefficient contains too many non-mesh dimensions
	error("ERROR: the coefficient contains too many dimensions that differ" + ...
		" in length from that of the mesh measures");
else
	% coefficient dimensions are incompatible with the allowed forms
	% error("ERROR: the coefficient has an invalid size") %UPDATE
end

% check rank
rank = numel(where_nDim);
if rank == 0
	d1_d2 = [ndims(data)+1, ndims(data)+2];
elseif rank == 1
	d1_d2 = [where_nDim, ndims(data)+1];
elseif rank ==2
	d1_d2 = where_nDim;
else
	error("ERROR: the coefficient appears to be a rank "+rank+...
		" tensor. Tensors with rank>2 are not supported")
end

if isConstant
	% pad the array dimensions if the field is constant
	where_othr = ndims(data)+(3-rank);
end

% permute the coefficient dimensions
dimOrder = [d1_d2, where_nPts, where_othr, where_nTri];
dimOrder = [dimOrder,setdiff(1:ndims(data),dimOrder)];
data = permute(data,dimOrder);
% DONE
end