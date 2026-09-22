function [S] = femat_surface(mesh,scalarC,vectorR,bndFaces)
% constructs the finite-element surface integral matrix
arguments
	mesh FE_mesh
	scalarC FE_coeff
	vectorR FE_coeff = grad(FE_coeff.baryfunc(mesh));
	bndFaces = mesh.boundaryFaces;
end

% check the coefficients
if (scalarC.rank == 1)
	error("ERROR: the second input must be a scalar/tensor field")
end
if (vectorR.rank ~= 1)
	error("ERROR: the second input must be a vector field")
elseif size(vectorR,2)~=1
	vectorR = vectorR.';
end

% get relevant portion of the triangulation
[numPts,numDim,numPoly,numTri] = mesh.num("pts","dim","poly","tri");
[~,numFaceVtx] = size(bndFaces);
faceNodes = unique(bndFaces(:));
triOnBnd = find(sum(ismember(mesh.tri,faceNodes),1)>=numFaceVtx);

% get the baryvectors of the nodes adjacent to each face
bvec = zeros(numDim,1,numPoly,numTri);
facePerms = ncombsk(numPoly,numFaceVtx,"asIndex",true);
bndFaces = sort(bndFaces,2);
for ii = 1:size(facePerms,1)
	perm_ii = facePerms(ii,:);
	faces_ii = sort(mesh.tri(perm_ii,triOnBnd),1).';
	hasFace_ii = triOnBnd(ismember(faces_ii,bndFaces,"rows"));

	bvec_ii = mean(mesh.baryvec(:,:,~perm_ii,hasFace_ii),3);
	bvec(:,:,perm_ii,hasFace_ii) = bvec(:,:,perm_ii,hasFace_ii) + bvec_ii;
end

% While the surface integral depends on the surface normals (nvec = bvec/|bvec|)
% and the surface/face areas, my integration algorithm is set up for integrals over 
% the mesh's full dimensionality. Thus we need to find apply a scaling factor to each
% normal-vector so that each simplex volumes is converted to a surface/face area. 
% Luckily, given an n-dimensional simplex, the area of the (n-1)-dimensional sub-simplex 
% is given by: (area = n * vol/h), where 'vol' is the simplex volume, and 'h' is the 
% height of the simplex in the direction normal to the (n-1) face. Given 'bvec',the
% baryvector of the vertex opposite to the (n-1) sub-simplex, the height 'h' is
% equal to (h = 1/|bvec|). For the conversion factor (vol2area = area/vol = n*|bvec|),
% the adjusted normals (vol2area * nvec) simplify to (n * bvec).

% adjusted normals
nvec = numDim * bvec;

% convert to coefficient
vectorL = FE_coeff(mesh,nvec,mesh.vtxPoly).';

% construct the matrix
if nargin==2
	% vectorR is the gradient of the baryfunctions
	[vals,ii,jj] = FE_coeff.integrate([1,3],vectorL,scalarC,vectorR);
	S = sparse(mesh.tri(ii,:),mesh.tri(jj,:),squeeze(vals),numPts,numPts);
else
	% vectorR is another coefficient
	bf = FE_coeff.baryfunc(mesh);
	[vals,ii,jj] = FE_coeff.integrate([1,4],vectorL,scalarC,vectorR,bf);
	S = sparse(mesh.tri(ii,:),mesh.tri(jj,:),squeeze(vals),numPts,numPts);
end
% DONE
end