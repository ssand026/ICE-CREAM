function [baryVec] = gbaryvectors(pts,tri)
% Returns the set of barycentric vectors for each simplex in the mesh
% The output is a [numDim x numVtx x numTri] array
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end
[numDim,~] = size(pts);
[numVtx,numTri] = size(tri);

% construct baryvectors
allPts = reshape(pts(:,tri),numDim,numVtx,numTri);
baryVec = zeros(numDim,numVtx,numTri);
for nn = 1:numVtx
	indx = [1:nn-1,nn+1:numVtx];

	% get side-length vectors
	sideVec = allPts(:,indx,:) - allPts(:,indx(end),:);
	sideVec(:,end,:) = 1;
	sideNorm = vecnorm(sideVec,2,1);
	sideVec = sideVec ./ (sideNorm + ~sideNorm);

	% get normal vectors (via Gram-Schmidt orthogonalization)
	Q = zeros(size(sideVec));
	for jj = 1:size(Q,2)
		Q(:,jj,:) = sideVec(:,jj,:);
		if jj > 1
			proj = sum(Q(:,1:jj-1,:).*Q(:,jj,:),1);
			proj = sum(Q(:,1:jj-1,:).*proj,2);
			Q(:,jj,:) = Q(:,jj,:) - proj;
		end
		% normalize
		nQ = vecnorm(Q(:,jj,:),2,1);
		Q(:,jj,:) = Q(:,jj,:) ./ (nQ + ~nQ);
	end
	normVec = Q(:,end,:);

	% project normals onto side length to get baryvectors
	lenVec = allPts(:,nn,:) - allPts(:,indx(1),:);
	baryVec(:,nn,:) = normVec ./ sum(normVec .* lenVec,1);
end
% DONE
end