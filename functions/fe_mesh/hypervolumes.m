function [vol] = hypervolumes(pts,tri)
% Returns the hypervolume of each simplex in the triangulation.
% Uses the Calyley-Menger determinant, so the triangulation can be of arbitrary
% dimensionality (computes edge lengths, triangle areas, tetrahedron volumes, etc.)
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end
numDim = size(pts,1);
numVtx = size(tri,1);
numTri = size(tri,2);

if numVtx<=1
	error("ERROR: the triangulation must contain two or more points per simplex")
elseif numVtx==2
	% simply compute the edge lengths
	vol = vecnorm(pts(:,tri(1,:))-pts(:,tri(2,:)),2,1);
else
	% construct the Cayley-Menger matrix
	S = zeros(numVtx+1,numVtx+1,numTri);
	S(end,:,:) = 1;
	S(:,end,:) = 1;
	S(end,end,:) = 0;
	
	vtx = reshape(pts(:,tri),numDim,numVtx,numTri);
	for ii = 1:numVtx
		for jj = (ii+1):numVtx
			distSqrd = sum((vtx(:,ii,:)-vtx(:,jj,:)).^2,1);
			S(ii,jj,:) = distSqrd;
			S(jj,ii,:) = distSqrd;
		end
	end
	% compute the determinant
	detT = reshape(prod(pageeig(S),1),1,numTri);
	
	% take the squre root and re-scale
	vol = sqrt(abs(detT)/(2^(numVtx-1))) / factorial(numVtx-1);
end
% DONE
end