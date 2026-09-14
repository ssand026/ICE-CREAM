function [vals] = expand(mesh,vals)
% expand -  Expands a coefficient so its value is defined on each simplex
%
% The resulting array will have dimensions:
%  - [d1, d2, 1, 1]:       if coeff is domain constant
%  - [d1, d2, 1, numTri]:    if coeff is simplex constant
%  - [d1, d2, numPoly, numTri]: if coeff is domain/simplex continuous
%
% SEE ALSO: FE_coeff
arguments
	mesh FE_mesh
	vals double
end

[vals,isDisjoint,isConstant] = check(mesh,vals);
if isDisjoint && ~isConstant
	% field is variable over each simplex
	% (do nothing)
	
elseif isDisjoint && isConstant
	% field is constant over each simplex
	[numPoly] = mesh.num("poly");
	vals = repmat(vals,1,1,numPoly,1);

elseif ~isDisjoint && ~isConstant
	% field is variable over the domain
	vals = nodal_to_vertex(mesh,vals);

elseif ~isDisjoint && isConstant
	% field is constant over the domain
	[numTri,numPoly] = mesh.num("tri","poly");
	vals = repmat(vals,1,1,numPoly,numTri);

end
% DONE
end