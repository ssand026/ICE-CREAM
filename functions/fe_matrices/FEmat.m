function [out] = FEmat(mesh, matrixType, scalarC, vectorA, vectorB)
% Returns the specified finite-element matrix
%
% INPUTS:
%	mesh: mesh for the system
%	matType: type of matrix to construct
%	scalarC: scalar/tensor coefficient for the matrix
%	vectorA: 1st vector coefficient for the matrix
%   vectorB: 2nd vector coefficient for the matrix
%
% SEE ALSO: FE_MESH, FE_COEFF
arguments
	mesh FE_mesh
	matrixType (1,1) {mustBeMember(matrixType,["overlap","stiffness","scalar", ...
		"skew-vector","symm-vector","dirichlet","neumann","surface"])}
	scalarC {mustBeA(scalarC,{'FE_coeff','function_handle','char','string','numeric'})} = [];
	vectorA {mustBeA(vectorA,{'FE_coeff','function_handle','char','string','numeric'})} = [];
	vectorB {mustBeA(vectorB,{'FE_coeff','function_handle','char','string','numeric'})} = [];
end

%% check the first coefficient
%===========================================================
if matches(matrixType,["scalar","stiffness","surface","skew-vector","symm-vector"])
	% attempt to evaluate the coefficient
	if ~isa(scalarC,"FE_coeff") && ~isempty(scalarC)
		try
			scalarC = FE_coeff(mesh,scalarC);
		catch ME
			errorMsg = "ERROR: the first coefficient does not evaluate to a valid FE_coeff";
			error(errorMsg + newline + ME.message);
		end
	end
	% check existance/rank
	if isempty(scalarC)
		error("ERROR: a scalar/tensor coefficient was not specified")
	elseif isa(scalarC,"FE_coeff") && (matrixType=="scalar") && (scalarC.rank~=0)
		error("ERROR: the coefficient must be a scalar quantity")
	elseif isa(scalarC,"FE_coeff") && (scalarC.rank==1)
		error("ERROR: the first coefficient must be a scalar/tensor quantity")
	end
end

%% check the second coefficient
%===========================================================
if matches(matrixType,["surface","skew-vector","symm-vector"])
	% attempt to evaluate the coefficient
	if ~isa(vectorA,"FE_coeff") && ~isempty(vectorA)
		try
			vectorA = FE_coeff(mesh,vectorA);
		catch ME
			errorMsg = "ERROR: the second coefficient does not evaluate to a valid FE_coeff";
			error(errorMsg + newline + ME.message);
		end
	end
	% check existance/rank
	if isempty(vectorA) && (matrixType~="surface")
		error("ERROR: the first vector coefficient was not specified")
	elseif isa(vectorA,"FE_coeff") && (vectorA.rank~=1)
		error("ERROR: the second coefficient must be a vector quantity")
	end
	
	% check the vector orientation
	if matches(matrixType,["surface","skew-vector"])
		% must be a column-vector
		if size(vectorA,1)==1; vectorA = vectorA.'; end
	else
		% must be a row-vector
		if size(vectorA,2)==1; vectorA = vectorA.'; end
	end
end

%% check the third coefficient
%===========================================================
if matches(matrixType,["symm-vector"])
	% attempt to evaluate the coefficient
	if ~isa(vectorB,"FE_coeff") && ~isempty(vectorB)
		try
			vectorB = FE_coeff(mesh,vectorB);
		catch ME
			errorMsg = "ERROR: the third coefficient does not evaluate to a valid FE_coeff";
			error(errorMsg + newline + ME.message);
		end
	end
	% check existance/rank
	if isempty(vectorB)
		error("ERROR: the second vector coefficient was not specified")
	elseif isa(vectorB,"FE_coeff") && (vectorB.rank~=1)
		error("ERROR: the third coefficient must be a vector quantity")
	end
	
	% check the vector orientation (must be a column-vector)
	if size(vectorB,1)==1; vectorB = vectorB.'; end
end

%% compute the matrix
%===========================================================
% define the test/basis-functions for any locally-assembled matrices
if matches(matrixType,["overlap","scalar","stiffness","skew-vector","symm-vector"])
	[numPts] = mesh.num("pts");
	tri = mesh.triPoly;
	bf = FE_coeff.baryfunc(mesh);
end

% select matrix type
switch matrixType
	case "dirichlet"
		% get the Dirichlet boundary-condition matrices
		[intoB,exitB] = femat_dirichlet(mesh);
		out = {intoB,exitB};
	case "neumann"
		% get the Neumann boundary-condition matrices
		[intoB,exitB] = femat_neumann(mesh);
		out = {intoB,exitB};
	case "orthobasis"
		% get the orthonormal dirichlet boundary-condition matrix
		[intoB,exitB] = femat_orthobasis(mesh);
		out = {intoB,exitB};
	case "surface"
		% get the surface integral matrix
		if isempty(vectorA)
			S = femat_surface(mesh,scalarC);
		else
			S = femat_surface(mesh,scalarC,vectorA);
		end
		out = S;
	case "overlap"
		% assemble the overlap/mass matrix
		[vals,ii,jj] = FE_coeff.integrate([1,2],bf,bf);
		M = sparse(tri(ii,:),tri(jj,:),squeeze(vals),numPts,numPts);
		M = (M + M.')/2;
		out = M;
	case "stiffness"
		% assemble the stiffness matrix
		df = grad(bf);
		[vals,ii,jj] = FE_coeff.integrate([1,3],df.',scalarC,df);
		K = sparse(tri(ii,:),tri(jj,:),squeeze(vals),numPts,numPts);
		K = (K + K.')/2;
		out = K;
	case "scalar"
		% assemble the scalar-potential matrix
		[vals,ii,jj] = FE_coeff.integrate([2,3],scalarC,bf,bf);
		V = sparse(tri(ii,:),tri(jj,:),squeeze(vals),numPts,numPts);
		V = (V + V.')/2;
		out = V;
	case "skew-vector"
		% assemble the skew-symmetric (1st order) vector-potential matrix
		df = grad(bf);
		[valsL,ii,~] = FE_coeff.integrate([1,2],bf,df.',scalarC,vectorA);
		[valsR,~,jj] = FE_coeff.integrate([1,2],df.',bf,scalarC,vectorA);
		vals = valsL - valsR;
		X = sparse(tri(ii,:),tri(jj,:),squeeze(vals),numPts,numPts);
		X = (X - X.')/2;
		out = X;
	case "symm-vector"
		% assemble the symmetric (2nd-order) vector-potential matrix
		[vals,ii,jj] = FE_coeff.integrate([4,5],vectorA,scalarC,vectorB,bf,bf);
		Z = sparse(tri(ii,:),tri(jj,:),squeeze(vals),numPts,numPts);
		Z = (Z + Z.')/2;
		out = Z;
end
% DONE
end