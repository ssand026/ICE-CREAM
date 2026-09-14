function [nodeData] = vertex_to_nodal(mesh,vtxData)
% converts a vertex-defined coefficient to a node-defined coefficient
% vertex (disjoint) coefficients explicitly specify the coefficient's
% value at each point in the simplex set, whereas node coefficients
% have a single value on each mesh node 
arguments
	mesh FE_mesh
	vtxData (:,:,:,:) double
end

% circumvent the costly overhead of "unique(...)" using:
%	[C,ia,ic] = matlab.internal.math.uniquehelper(A,doSort,isFirst,byRows)
fast_unique = @(x,opts) matlab.internal.math.uniquehelper(x,opts(1),opts(2),opts(3));

% get the polynomial interpolation points
[numPts,numTri,numPoly] = mesh.num("pts","tri","poly");
tri = mesh.triPoly;
tri = double(tri(:));

% convert singular vertex data
if isequal(size(vtxData),[numPoly,numTri])
	node_value = fast_unique([tri,vtxData(:)],[false,true,true]);
	count = accumarray(node_value(:,1),1,[numPts,1],[],1);
	total = accumarray(node_value(:,1),node_value(:,2),[numPts,1]);
	nodeData = total./count;
	return
end

% convert FE_coeff style data
[numRow,numCol] = size(vtxData,[1,2]);
num_ij = (numRow*numCol);
[ii,jj] = ind2sub([numRow,numCol],1:num_ij);

nodeData = zeros(numRow,numCol,numPts);
for nn = 1:num_ij
	vtx_ij = vtxData(ii(nn),jj(nn),:);
	node_value = fast_unique([tri,vtx_ij(:)],[false,true,true]);
	count = accumarray(node_value(:,1),1,[numPts,1],[],1);
	total = accumarray(node_value(:,1),node_value(:,2),[numPts,1]);
	nodeData(ii(nn),jj(nn),:) = total./count;
end
% DONE
end