function [outerC,outerR] = outer_sphere(pts)
% returns the center/radius of the smallest sphere that encloses a set of points
arguments
	pts (:,:) double
end
pts = pts.';

% get points on the convex hull
[hullFaces,~] = convhulln(pts);
onHull = unique(hullFaces(:));
pts = pts(onHull,:);

% check if the sphere enclosing the two furthest points encloses all points
%-----------------------------------------------------------
[numPts,numDim] = size(pts);
sqrdDist = zeros(numPts,numPts);
for ii = 1:numDim
	sqrdDist = sqrdDist + (pts(:,ii).' - pts(:,ii)).^2;
end
[sqrdDiam,index] = max(sqrdDist,[],"all");
outerR = sqrt(sqrdDiam)/2;
[aa,bb] = ind2sub(size(sqrdDist),index);
outerC = (pts(aa,:) + pts(bb,:))/2;

if all(sum((pts - outerC).^2,2) <= (outerR+eps(outerR))^2)
	% we have found the circumscribing circle
	return
end

% the center of the enclosing sphere should lie on a vertex of the furthest-point voroni diagram
%-----------------------------------------------------------
% get the vertices of the furthest-point voronoi diagram
[vtx,~] = voronoin(pts,{'Qu','Tv'});
vtx = vtx(all(isfinite(vtx),2),:);

% find the largest distance between the point set and each vertex
maxDist = zeros(size(vtx,1),1);
for ii = 1:size(vtx,1)
	maxDist(ii,:) = sqrt(max(sum((vtx(ii,:) - pts).^2,2)));
end

% return the vertex with the smallest distance from its furthest neighbor
[outerR,index] = min(maxDist);
outerC = vtx(index,:);
end