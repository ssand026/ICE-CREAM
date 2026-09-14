function [vtxData] = nodal_to_vertex(mesh,nodeData)
% converts a node-defined coefficient to a vertex-defined coefficient
% vertex (disjoint) coefficients explicitly specify the coefficient's
% value at each point in the simplex set, whereas node coefficients
% have a single value on each mesh node 
arguments
	mesh FE_mesh
	nodeData (:,:,:) double
end

% get the polynomial interpolation points
[numPts,numTri,numPoly] = mesh.num("pts","tri","poly");
tri = mesh.triPoly;

% convert singular node-data
if isvector(nodeData) && (numel(nodeData)==numPts)
	vtxData = nodeData(tri);
	return
end

% convert FE_coeff style data
[numRow,numCol] = size(nodeData,[1,2]);
vtxData = zeros(numRow,numCol,numPoly,numTri);
for ii = 1:numRow
	for jj = 1:numCol
		node_ij = nodeData(ii,jj,:);
		node_ij = node_ij(:);
		vtxData(ii,jj,:,:) = node_ij(tri);
	end
end
% DONE
end