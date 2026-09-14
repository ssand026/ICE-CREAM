function [ptsPoly,triPoly,vtxConn,vtxPoly] = polymesh_params(pts,tri,degree)
% Returns the metrics needed to describe meshes containing polynomial basis functions
arguments
	pts (:,:)
	tri (:,:)
	degree (1,1) {mustBeInteger}
end
[~,numPts] = size(pts);
[numVtx,numTri] = size(tri);

% define the degree-1 components
ptsPoly = pts;
triPoly = tri;
vtxConn = speye(numPts);
vtxPoly = speye(numVtx);

for nn = 2:degree
	% get all unique faces of degree nn
	comb_nn = ncombsk(numVtx,nn,"asIndex",false);
	faces = reshape(tri(comb_nn.',:),nn,[]).';
	faces = unique(faces,"rows");
	numFace = size(faces,1);

	% get the point connectivity for the new nodes
	col = cast(repmat(1:numFace,nn,1).',"like",faces);
	vtxConn_nn = sparse(faces,col,true,numPts,numFace);

	% get the barycentric powers for the new vertices
	vtxPoly_nn = ncombsk(numVtx,nn,"asIndex",true);

	% get the location of the polynomial points
	ptsPoly_nn = full(pts * vtxConn_nn)./full(sum(vtxConn_nn,1));

	% get the triangulations for the polynomial points
	numCombs = size(comb_nn,1);
	triPoly_nn = zeros(numCombs,numTri,"like",tri);
	for ii = 1:numCombs
		face_ii = tri(comb_nn(ii,:),:).';
		[~,triPoly_nn(ii,:)] = ismember(face_ii,faces,"rows");
	end
	triPoly_nn = triPoly_nn + size(ptsPoly,2);

	% update polynomial mesh properties
	ptsPoly = [ptsPoly, ptsPoly_nn]; %#ok
	triPoly = [triPoly; triPoly_nn]; %#ok
	vtxConn = [vtxConn, vtxConn_nn]; %#ok
	vtxPoly = [vtxPoly; vtxPoly_nn]; %#ok
end
% DONE
end
