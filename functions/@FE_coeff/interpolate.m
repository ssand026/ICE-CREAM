function [values] = interpolate(obj,queryPts)
% returns the value of the coefficient at the given set of points
arguments
	obj FE_coeff
	queryPts (:,:) double
end

% check input validity
if size(queryPts,1) ~= obj.mesh.num("dim")
	error("ERROR: the given set of points must have "+obj.mesh.num("dim")+" rows.")
end


if obj.isConstant && ~obj.isDisjoint
	% coefficient is constant over the domain
	values = repmat(obj.vals,1,1,size(queryPts,2));
	values = squeeze(values);
	return
end

% get the enclosing triangle and barycentric coordinates of each query point
TR = triangulation(double(obj.mesh.tri.'), obj.mesh.pts.');
[queryTri,barycoord] = TR.pointLocation(queryPts.');

if obj.isConstant && obj.isDisjoint
	% coefficient is constant over each simplex
	values = obj.vals(:,:,:,queryTri);
	values = squeeze(values);
	return
end

% get the polynomial interpolation weights
polyWeights = full(obj.poly) * barycoord.';

% multiply the interpolation weights by the coefficient polynomial terms
% to get the real value of the coefficient at the given points
values = tensormult(polyWeights,'ct',obj.vals(:,:,:,queryTri),'abct','abt');
values = squeeze(values);
% DONE
end