function [coeffs,indx,poly] = baryintegral(vtxPoly,numBary)
% Compute the indices and integration coefficients for the product of 'n'
% barycentric test-functions, where each test-functions is a polynomial defined
% by the sum of barycoordinate powers in each row of the array 'vtxPoly'.
arguments
	vtxPoly (:,:) double {mustBeGreaterThan(vtxPoly,-1),mustBeFinite}
	numBary (1,1) double {mustBeInteger,mustBePositive}
end

% get the polynomial powers for a single barycentric basis-function
numBasis = size(vtxPoly,1);
basisPoly = vtxPoly;
basisIndx = (1:numBasis).';

% compute the indexed product of the basis-functions
totalPoly = basisPoly;
totalIndx = basisIndx;
for ii = 2:numBary
	numTotal = size(totalPoly,1);
	totalPoly = repmat(totalPoly,numBasis,1) + repelem(basisPoly,numTotal,1);
	totalIndx = [repmat(totalIndx,numBasis,1), repelem(basisIndx,numTotal,1)];
end

% Compute the coefficients of the barycentric integral from the polynomial powers.
% Multiply by the simplex volumes to obtain the full integration coefficients.
poly = totalPoly;
indx = totalIndx;
numDim = size(poly,2) - 1;
coeffs = factorial(numDim) * (prod(gamma(1+poly),2) ./ gamma(1+numDim+sum(poly,2)));

% DONE
end
