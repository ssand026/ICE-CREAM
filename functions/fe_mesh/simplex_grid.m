function [baryc,tri] = simplex_grid(numVtx,numSeg)
% returns the barycoordinates and triangulations for the set of 
% points obtained from regular sub-division of a simplex region 
arguments
	numVtx double {mustBeInteger,mustBeGreaterThan(numVtx,1)}
	numSeg double {mustBeInteger,mustBePositive}
end

if (numSeg == 1)
	% no subdivision of the simplex, return the base triangulation
	baryc = eye(numVtx);
	tri = (1:numVtx).';
	return
end

% generate the barycoordinates for a simplex where each side has been 
% split into 'numSeg' segments of equal length
baryc = (0:numSeg).';
for ii = 2:numVtx
	baryc = [repmat(baryc,numSeg+1,1), repelem((0:numSeg).',size(baryc,1),1)];
	baryc((sum(baryc,2) > numSeg),:) = [];
end
baryc = baryc((sum(baryc,2) == numSeg),:);

% convert to cartesian coordinates
coords = baryc * [zeros(1,numVtx-1); eye(numVtx-1)];

% normalize barycoordinates
baryc = baryc / numSeg;

if (nargout < 2)
	% exit without calculating the triangulation
	tri = [];
	return
end

% get distance between node pairs
numGrid = size(coords,1);
nDist = zeros(numGrid,numGrid);
for dd = 1:numVtx-1
	nDist = nDist + abs(coords(:,dd)-coords(:,dd).');
end

% get absolute distance from the origin
aDist = sum(coords,2);

% get the adjacency matrix
isAdjacent = triu(nDist==1);
isDiagonal = ((nDist == 2) & (aDist(:)==aDist(:).'));
adjMat = (isAdjacent | isDiagonal);
adjMat = (adjMat | adjMat.');

% create the edge list
[row,col] = find(triu(adjMat));
edges = unique([row,col],"rows");

% create triangulation from the adjacency matrix
G = graph(edges(:,1),edges(:,2));
tri = allcycles(G,'MinCycleLength',numVtx,'MaxCycleLength',numVtx);
tri = reshape([tri{:}].',numVtx,[]);
tri = unique(sort(tri,1).',"rows").';

% remove cycles containing invalid edges
edgePerms = ncombsk(numVtx,2);
hasEdge = true(size(tri,2),1);
for ii = 1:size(edgePerms,1)
	edges_ii = tri(edgePerms(ii,:),:);
	hasEdge = hasEdge & ismember(edges_ii.',edges,"rows");
end
tri = tri(:,hasEdge);
% DONE
end