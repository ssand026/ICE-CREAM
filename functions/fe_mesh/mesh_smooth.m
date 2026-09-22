function [pts] = mesh_smooth(pts,tri,reg,opt)
% Improves the quality of a mesh by shifting nodes towards the center of their neighbors
% The input "iters" determines how many times the node positions are re-adjusted,
% with larger values resulting in more smoothing of the mesh.
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	reg (1,:) {mustBeInteger,mustBeNonnegative} = ones(1,size(tri,2));
	opt.iters (1,1) {mustBeInteger,mustBePositive} = 10;
	opt.fixedSurfs (1,1) {mustBeMember(opt.fixedSurfs, ...
		["all","domain","region","none"])} = "all"; % set of surfaces that are allowed to move
	opt.tolerance (1,1) {mustBeInRange(opt.tolerance,0,1)} = 0.99; % tolerance for orthoganality condition
end

[numDim,numPts] = size(pts);

% find faces/nodes on each boundary surface
[bndFaces,bndRegions] = region_bounds(tri,reg);
bndNodes = cellfun(@(x) unique(x(:)).',bndFaces,"UniformOutput",false);

% find the set of nodes unique to each boundary surface
[allBndNodes,~,ic] = unique([bndNodes{:}]);
multiBndNodes = allBndNodes(accumarray(ic(:),1)>1);
uniqueNodes = cellfun(@(x) setdiff(x,multiBndNodes),bndNodes,"UniformOutput",false);

% find any nodes without constraints on their position
%-----------------------------------------------------------
canMoveFree = true(1,numPts);
canMoveFree(allBndNodes) = false;

% find any nodes that can only move orthogonally to the boundaries
%-----------------------------------------------------------
canMoveOrth = false(1,numPts);
normVecs = zeros(numDim,numPts);

% get the set of free boundary elements
switch opt.fixedSurfs
	case "all"
		freeSurfs = [];
	case "region"
		freeSurfs = find(cellfun(@(x) any(x==0), bndRegions));
	case "domain"
		freeSurfs = find(cellfun(@(x) ~any(x==0), bndRegions));
	case "none"
		freeSurfs = [1:numel(bndFaces)];
end

for ii = freeSurfs(:).'
	% find the normals for faces on the boundary
	if isempty(reg) || isscalar(reg)
		[nvec_ii] = boundary_normals(pts,tri,bndFaces{ii});
	else
		bndReg_ii = bndRegions{ii};
		bndFace_ii = bndFaces{ii};
		tri_ii = tri(:,reg==max(bndReg_ii));
		[nvec_ii] = boundary_normals(pts,tri_ii,bndFace_ii);
	end

	% find nodes whose normals all agree within tol
	for jj = uniqueNodes{ii}
		% get the normals of faces containing node jj
		face_jj = any(bndFaces{ii}==jj,2);
		nvec_jj = nvec_ii(:,face_jj);
		nvec_avg = mean(nvec_jj,2);

		% check if normal vectors agree within tolerance
		overlap = abs(nvec_avg.' * nvec_ii);
		if all(overlap >= opt.tolerance,"all")
			normVecs(:,jj) = nvec_avg;
			canMoveOrth(jj) = true;
		end
	end
end
% remove any unused normals
normVecs = normVecs(:,canMoveOrth);

% begin shifting the node locations
%-----------------------------------------------------------
if ~any(canMoveFree) && ~any(canMoveOrth)
	% all nodes are fixed in place, exit before loop
	return
end

% normalize and pad the adjacency matrix
avgMat = adjacency_matrix(tri);
avgMat = avgMat ./ (~sum(avgMat,1) + sum(avgMat,1));
avgMat(end+1:numPts,:) = 0;
avgMat(:,end+1:numPts) = 0;

for nn = 1:opt.iters
	% for each node, get the average position its neighbors
	newPos = (pts * avgMat);

	% update any un-constrained nodes
	if any(canMoveFree)
		pts(:,canMoveFree) = newPos(:,canMoveFree);
	end

	% update any othogonally-constrained nodes
	if any(canMoveOrth)
		shiftDelta = pts(:,canMoveOrth) - newPos(:,canMoveOrth);
		projNorm = dot(normVecs,shiftDelta).*normVecs;
		pts(:,canMoveOrth) = pts(:,canMoveOrth) - (shiftDelta - projNorm);
	end
end
% DONE
end