function [innerC,innerR] = inner_sphere(pts,TR)
% returns the center/radius of the largest sphere enclosed by a set of points
arguments
	pts (:,:) double {mustBeReal}
	TR {mustBeA(TR,{'triangulation','delaunayTriangulation'})} = delaunayTriangulation(pts.');
end
pts = pts.';

% get the vertices of the nearest-point voronoi diagram
[vtx,~] = voronoin(pts);
vtx = vtx(all(isfinite(vtx),2),:);
vtx = vtx(~isnan(pointLocation(TR,vtx)),:);
vtx = vtx(~ismember(vtx,pts,"rows"),:);

% find the smallest distance between the point set and each vertex
minDist = zeros(size(vtx,1),1);
for ii = 1:size(vtx,1)
	minDist(ii,:) = sqrt(min(sum((vtx(ii,:) - pts).^2,2)));
end

% return the vertex with the largest distance from its' closest neighbor
[innerR,largestInner] = max(minDist);
innerC = vtx(largestInner,:);
end