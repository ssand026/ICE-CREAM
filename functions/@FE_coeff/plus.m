function [C] = plus(A,B)
% computes the sum of two scalar/vector/tensor fields
arguments
	A {mustBeA(A,{'FE_coeff','numeric','logical'})}
	B {mustBeA(B,{'FE_coeff','numeric','logical'})}
end

isCoeffA = isa(A,"FE_coeff");
isCoeffB = isa(B,"FE_coeff");

% compute the sum of ...
if (isCoeffA && isCoeffB)
	% check if coeffs use the same mesh
	if A.mesh ~= B.mesh
		error("ERROR: the FE_coeffs are defined on different meshes")
	else
		meshC = A.mesh;
	end
elseif (isCoeffA && ~isCoeffB)
	% convert B to coeff
	try
		meshC = A.mesh;
		B = FE_coeff(meshC,B);
	catch ME
		error("ERROR: coefficients have incompatible sizes." + newline + ME.message)
	end
elseif (~isCoeffA && isCoeffB)
	% convert A to coeff
	try
		meshC = B.mesh;
		A = FE_coeff(meshC,A);
	catch ME
		error("ERROR: coefficients have incompatible sizes." + newline + ME.message)
	end
else
	% it shouldn't be possible to end up here
	error("ERROR: at least one input must be a FE_coeff object.")
end

% circumvent the costly overhead of "unique(...)" using:
%	[C,ia,ic] = matlab.internal.math.uniquehelper(A,doSort,isFirst,byRows)
fast_unique = @(x,opts) matlab.internal.math.uniquehelper(x,opts(1),opts(2),opts(3));

% find identical polynomial terms
[polyC,~,CtoAB] = fast_unique(full([A.poly; B.poly]),[false,true,true]);

% initialize the combined data
sizeC = max(size(A.vals,1:4), size(B.vals,1:4));
dataC = zeros(sizeC(1),sizeC(2),size(polyC,1),sizeC(4));

lenA = size(A.poly,1);
lenB = size(B.poly,1);

% add contributions from coeff A
for ii = 1:lenA
	c_ii = CtoAB(ii);
	dataC(:,:,c_ii,:) = dataC(:,:,c_ii,:) + A.vals(:,:,ii,:);
end

% add contributions from coeff B
for ii = 1:lenB
	c_ii = CtoAB(lenA+ii);
	dataC(:,:,c_ii,:) = dataC(:,:,c_ii,:) + B.vals(:,:,ii,:);
end

% construct new coeff
C = FE_coeff(meshC, dataC, polyC);

% DONE
end