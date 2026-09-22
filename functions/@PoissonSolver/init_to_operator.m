function [scalar2operator] = init_to_operator(mesh)
% INIT_TO_OPERATOR Returns a function that converts scalar 
% coefficients into their operator (matrix) form
arguments
	mesh FE_mesh
end

% load common parameters
[numPts,numTri,numPoly] = mesh.num("pts","tri","poly");
vol = mesh.volumes;
tri = mesh.triPoly;

% compute integral coefficients
[int_ijk] = baryintegral(mesh.vtxPoly,3);
int_ijk = reshape(int_ijk,numPoly,numPoly,numPoly);

% accumulate over row/column indices
[ii,jj] = ind2sub([numPoly,numPoly],[1:numPoly^2]);
t_ii = reshape(tri(ii,:),[],1);
t_jj = reshape(tri(jj,:),[],1);

[tri_ij,~,accumIndx] = unique([t_ii,t_jj],"rows");
t_ii = tri_ij(:,1);
t_jj = tri_ij(:,2);
numUnique = size(tri_ij,1);

% return the converter function
scalar2operator = @scalar_to_operator;

% internal functions
%***********************************************************
	function [V] = scalar_to_operator(coeff)
	% reshape the input coefficient
	coeff = squeeze(coeff);
	sz = size(coeff);
	% apply the simplex volumes
	if isequal(sz,[numPts,1]) || isequal(sz,[1,numPts])
		tmp = vol.*coeff(tri);
	elseif isequal(sz,[numPoly,numTri])
		tmp = vol.*(coeff);
	elseif isequal(sz,[numTri,numPoly])
		tmp = vol.*(coeff.');
	else
		error("ERROR: input coefficient has invalid size")
	end
	% convert to operator
	tmp = tensormult(int_ijk,'ijk',tmp,'kt','ijt'); % apply integral coeffs
	tmp = accumarray(accumIndx,tmp(:),[numUnique,1]); % accumulate identical subscripts
	V = sparse(t_ii,t_jj,tmp,numPts,numPts);
	V = (V + V.')/2;
	end
%***********************************************************
end