function [wavefunc2rhs,density2rhs,braket2rhs] = init_to_rhs(mesh)
% INIT_TO_RHS returns the functions that construct the right side of
% the Poisson equation from a set of charge-densities, wavefunctions, or
% bra-ket pairs.
arguments
	mesh FE_mesh
end

% load common parameters
[numPts,numPoly,numTri] = mesh.num("pts","poly","tri");
vol = mesh.volumes;
tri = mesh.triPoly;

% compute the integration coefficients if inputs are charge-densities 
% (rhs is first-order on the domain)
[int_ab] = baryintegral(mesh.vtxPoly,2);
int_ab = reshape(int_ab,numPoly,numPoly);

% compute the integration coefficients if inputs are wavefunctions
% (rhs is second-order on the domain)
[int_ijk] = baryintegral(mesh.vtxPoly,3);
int_ijk = reshape(int_ijk,numPoly,numPoly,numPoly);

% return the converter functions
wavefunc2rhs = @wavefunc_to_rhs;
density2rhs = @density_to_rhs;
braket2rhs = @braket_to_rhs;

% internal functions
%***********************************************************
	function [rhs] = density_to_rhs(rho,dim)
	% computes the right side of the Poisson equation from a charge-density

	% reshape input
	rho = data_to_vertex(rho,tri,[numPts,numPoly,numTri],dim);
	numRho = size(rho,3);
	
	% return the rhs
	tmp = tensormult(int_ab,'ab',rho,'bts','ats');
	tmp = reshape(vol .* tmp,[],numRho);
	rhs = zeros(numPts,numRho);
	for ss = 1:numRho
		rhs(:,ss) = accumarray(tri(:),tmp(:,ss),[numPts,1]);
	end
	% DONE
	end
%***********************************************************
	function [rhs] = wavefunc_to_rhs(psi,dim)
	% computes the right side of the Poisson equation from a wavefunction
	
	% reshape input
	ket = data_to_vertex(psi,tri,[numPts,numPoly,numTri],dim);
	numKet = size(ket,3);

	% return the rhs
	tmp = tensormult(int_ijk,'ijk',ket,'jts','ikts');
	tmp = tensormult(conj(ket),'its',tmp,'ikts','kts');
	tmp = reshape(vol .* tmp,[],numKet);
	rhs = zeros(numPts,numKet);
	for ss = 1:numKet
		rhs(:,ss) = accumarray(tri(:),tmp(:,ss),[numPts,1]);
	end
	% DONE
	end
%***********************************************************
	function [rhs] = braket_to_rhs(bra,ket,dim)
	% computes the right side of the Poisson equation from a pair of bra/kets
	
	% check the 'dim' input
	if isempty(dim)
		% assume bra/ket have a set size of one
		braDim = [];
		ketDim = [];
	elseif isscalar(dim)
		% assume sets lie along the same dimension in bra/ket
		braDim = dim;
		ketDim = dim;
	elseif isequal(size(dim),[1,2])
		% sets dimensions for bra/ket were given
		braDim = dim(1);
		ketDim = dim(2);
	else
		% dim must have two or fewer elements
		error("ERROR: the input 'dim' cannot contain more than two elements")
	end

	% reshape bra/ket
	bra = data_to_vertex(bra,tri,[numPts,numPoly,numTri],braDim);
	ket = data_to_vertex(ket,tri,[numPts,numPoly,numTri],ketDim);
	numBra = size(bra,3);
	numKet = size(ket,3);

	% check that set dimensions are valid
	if ~isequal(numBra,numKet) && (numBra~=1) && (numKet~=1)
		error("ERROR: the 'bra'/'ket' inputs have incompatible set lengths")
	end

	% return the rhs
	numSets = max(numBra,numKet);
	tmp = tensormult(int_ijk,'ijk',ket,'jts','ikts');
	tmp = tensormult(bra,'its',tmp,'ikts','kts');
	tmp = reshape(vol .* tmp,[],numSets);
	rhs = zeros(numPts,numSets);
	for ss = 1:numSets
		rhs(:,ss) = accumarray(tri(:),tmp(:,ss),[numPts,1]);
	end
	% DONE
	end
%***********************************************************
% DONE
end

% external functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = data_to_vertex(data,tri,sizes,dim)
% restructures the data to be of size [numVtx,numTri,numData]
numPts = sizes(1);
numVtx = sizes(2);
numTri = sizes(3);

if isempty(dim)
	% assume numData = 1
	data = squeeze(data);
	sz = size(data);
	if isequal(sz,[numPts,1]) || isequal(sz,[1,numPts])
		data = data(tri);
	elseif isequal(sz,[numVtx,numTri])
		% data = data;
	elseif isequal(sz,[numTri,numVtx])
		data = data.';
	else
		error("ERROR: the input 'data' has invalid size")
	end
else
	numData = size(data,dim);
	basisDim = [1:dim-1,dim+1:ndims(data)];
	basisLen = size(data,basisDim);
	
	if isequal(basisLen,numPts)
		% data is node-data
		if isequal(basisDim,1)
			data = reshape(data(tri,:),numVtx,numTri,numData);
		elseif isequal(basisDim,2)
			data = reshape(data(:,tri),numData,numVtx,numTri);
			data = permute(data,[2,3,1]);
		else
			error("ERROR: the specified dimension is invalid")
		end
	elseif isequal(basisLen,[numVtx,numTri]) || isequal(basisLen,[numTri,numVtx])
		% data is vertex-data
		vtxDim = basisDim(basisLen==numVtx);
		triDim = basisDim(basisLen==numTri);
		data = permute(data,[vtxDim,triDim,dim]);
	else
		error("ERROR: the input has invalid size")
	end
end
% DONE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%