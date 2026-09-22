function [intoB,exitB] = femat_neumann(mesh,faces)
% Finite element basis reduction/expansion matrices for Neumann boundary-conditions
%
% The output 'intoB' is the basis reduction matrix
% The output 'exitB' is the basis expansion matrix
%
% Calculations done in the reduced basis will automatically satisfy Neumann
% boundary conditions on the given faces.
arguments
	mesh FE_mesh
	faces (:,:) {mustBeInteger,mustBePositive} = mesh.boundaryFaces;
end
[numDim,numPts,numVtx] = mesh.num("dim","pts","vtx");

% check the set of boundary faces
faces = unique(sort(faces,2),"rows");
if ~(size(faces,1) >= 1)
	error("ERROR: the inputted set of boundary faces is empty")
elseif (size(faces,2) ~= numDim)
	error("ERROR: the inputted set of boundary faces has invalid size")
end
numFaces = size(faces,1);

% extract the relevant portion of the triangulation
bndNodes = unique(faces(:));
isBndTri = (sum(ismember(mesh.tri,bndNodes),1)>=numDim);
tri = mesh.tri(:,isBndTri);
bvec = mesh.baryvec(:,:,:,isBndTri);

% --- construct the neumann-boundary matrix ---
N = sparse(numFaces,numPts);

% iterate over each permutation of the triangulation indices
facePerms = ncombsk(numVtx,numDim,"asIndex",true);
numPerms = size(facePerms,1);

lastRow = 0;
for ii = 1:numPerms
	% find faces on the boundary surface
	perm_ii = facePerms(ii,:);
	face_ii = sort(tri(perm_ii,:),1);
	onSurf = ismember(face_ii.',faces,"rows").';

	% find the surface normals
	adjVec = bvec(:,:,~perm_ii,onSurf);
	surfNorm = adjVec ./ sum(adjVec.^2,1);

	% add to boundary matrix
	val = squeeze(sum(surfNorm .* bvec(:,:,:,onSurf),1));
	rows = repmat(lastRow + (1:sum(onSurf)),numVtx,1);
	cols = tri(:,onSurf);
	N = N + sparse(cast(rows,"like",cols),cols,val,numFaces,numPts);
	
	% update filled row index
	lastRow = rows(end);
end

% if 'k' and 'd' are the indices of basis functions to keep/drop after applying
% the boundary conditions, we require any full-basis vector 'u', to meet the
% boundary conditions given by (N * u = g). Separating 'k' and 'd' yields:
%   [N(:,k), N(:,d)] * [u(k); u(d)] = g
% Rearranged
%   u(k) = N(:,k) \ (g - N(:,d) * u(d))
%   u(d) = N(:,d) \ (g - N(:,k) * u(k))
% Here we want to find the matrices 'intoB' and 'exitB' such that:
%   (intoB * u(k;d)) = u(k),  
%   (exitB * u(k)) = u(k;d),
% if g==0:
%   intoB = [eye(nk,nk),zeros(nk,nd)]
%   exitB = [eye(nk,nk), -N(:,d)\N(:,k)]
dropIndex = find(sum(mesh.vtxConn(bndNodes,:),1) == sum(mesh.vtxConn,1)); % bndNodes was nodes
keepIndex = setdiff((1:numPts).',dropIndex);

% get the Neumann boundary-condition application matrix
intoB = sparse(1:numel(keepIndex),keepIndex,1,numel(keepIndex),numPts);

% get the Neumann boundary-condition reversal matrix
[~,revertIndex] = sort([keepIndex; dropIndex]);
exitB = [speye(numel(keepIndex)); -(N(:,dropIndex) \ N(:,keepIndex))];
exitB = exitB(revertIndex,:);
% DONE
end
