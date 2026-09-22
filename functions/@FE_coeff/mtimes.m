function [C] = mtimes(A,B)
% computes the product between scalar/vector/tensor fields and/or a numeric array
arguments
	A {mustBeA(A,{'FE_coeff','numeric','logical'})}
	B {mustBeA(B,{'FE_coeff','numeric','logical'})}
end
isCoeffA = isa(A,"FE_coeff");
isCoeffB = isa(B,"FE_coeff");

% compute the product of ...
if (isCoeffA && isCoeffB)
	% coeff(A) * coeff(B)
	%-----------------------------------------------------------
	% check if coeffs use the same mesh
	if A.mesh ~= B.mesh
		error("ERROR: the FE_coeffs are defined on different meshes")
	else
		meshC = A.mesh;
	end
	
	% compute the tensor product
	if all(size(A.vals,[1,2])==1) || all(size(B.vals,[1,2])==1)
		% scalar products
		indicesA = 'ijac';
		indicesB = 'ijbc';
	else
		% vector/matrix products
		indicesA = 'ikac';
		indicesB = 'kjbc';
	end
	dataAB = tensormult(A.vals,indicesA,B.vals,indicesB,'ijabc');
	dataAB = reshape(dataAB,size(dataAB,1),size(dataAB,2),[],size(dataAB,5));
	
	% compute the polynomial-term products
	numA = size(A.poly,1);
	numB = size(B.poly,1);
	polyAB = repmat(A.poly,numB,1) + repelem(B.poly,numA,1);
	
	% check if either coefficient has stored indices
	if ~isempty(A.index) && ~isempty(B.index)
		indxAB = [repmat(A.index,numB,1), repelem(B.index,numA,1)];
	elseif ~isempty(A.index)
		indxAB = repmat(A.index,numB,1);
	elseif ~isempty(B.index)
		indxAB = repelem(B.index,numA,1);
	else
		indxAB = [];
	end
	
	% collect identical terms
	[dataC,polyC,indxC] = collect(dataAB,polyAB,indxAB);
	
	% construct new coeff
	C = FE_coeff(meshC, dataC, polyC, indxC);
	
elseif (isCoeffA && ~isCoeffB)
	% coeff(A) * array(B)
	%-----------------------------------------------------------
	if isscalar(A) || isscalar(B)
		dataC = (A.vals) * B;
	else
		dataC = tensormult(A.vals,'ijab',B,'jkab','ikab');
	end
	C = FE_coeff(A.mesh, dataC, A.poly);
	
elseif (~isCoeffA && isCoeffB)
	% array(A) * coeff(B)
	%-----------------------------------------------------------
	if isscalar(A) || isscalar(B)
		dataC = A * (B.vals);
	else
		dataC = tensormult(A,'ijab',B.vals,'jkab','ikab');
	end
	C = FE_coeff(B.mesh, dataC, B.poly);
	
else
	% array(A) * array(B) (it shouldn't be possible to end up here)
	%-----------------------------------------------------------
	error("ERROR: at least one input must be a FE_coeff object.")
end
% DONE
end