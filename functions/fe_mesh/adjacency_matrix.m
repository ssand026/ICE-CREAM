function [adjM] = adjacency_matrix(tri)
% Constructs the node adjacency matrix for the given triangulation.
arguments
	tri (:,:) {mustBeInteger,mustBePositive}
end
numVtx = size(tri,1);
numPts = max(tri,[],"all");

% estimate the number of nonzero elements
numNonzero = (numVtx^2-numVtx)*numPts;

% construct the adjacency matrix
adjM = logical(spalloc(numPts,numPts,numNonzero));
[ptA,ptB] = ncombsk(numVtx,2);
for nn = 1:numel(ptA|ptB)
	adjM = adjM | sparse(tri(ptA(nn),:),tri(ptB(nn),:),true,numPts,numPts);
end
adjM = (adjM.' | adjM);
% DONE
end