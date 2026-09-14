function [faces] = subfaces(tri,degree)
% Returns a list of all faces of the specified degree contained within the triangulation.
%
% If degree equals ...
%	1: returns the set of nodes in the mesh
%	2: returns the set of edges in the mesh
%	3: returns the set of triangles in the mesh
%	4: returns the set of tetrahedra in the mesh
%	etc.
%
arguments
	tri (:,:) {mustBeInteger,mustBePositive}
	degree (1,1) {mustBePositive,mustBeInteger}
end
[numVtx,~] = size(tri);

% circumvent the costly overhead of "unique(...)" using:
%	[C,ia,ic] = matlab.internal.math.uniquehelper(A,doSort,isFirst,byRows)
fast_unique = @(x,opts) matlab.internal.math.uniquehelper(x,opts(1),opts(2),opts(3));

if (degree > numVtx)
	% no sub-simplexes exist
	warning("The specified degree is larger than the triangulation size")
	faces = [];
elseif degree == 1
	% the set of all vertices in the triangulation
	faces = fast_unique(tri(:),[false,true,true]);
else
	% general solution
	tri = sort(tri,1);
	facePerms = ncombsk(numVtx,degree);
	faces = reshape(tri(facePerms.',:),degree,[]).';
	faces = fast_unique(faces,[false,true,true]);
end
% DONE
end