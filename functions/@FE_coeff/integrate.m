function [vals,varargout] = integrate(indexed,coeffs)
% INTEGRATE - integrate the product of several FE_coeffs
arguments
	indexed (1,:) double {mustBeInteger,mustBePositive}
end
arguments (Repeating)
	coeffs FE_coeff
end

% check that all coeffs are defined on the same mesh
meshes = cellfun(@(x)x.mesh,coeffs,"UniformOutput",false);
if ~isequal(meshes{:})
	error("ERROR: coefficients are defined on different meshes")
end

% check the validity of the kept index list
if max(indexed) > numel(coeffs)
	error("ERROR: the list of non-summed indices exceeds the number of coefficient terms ")
end

% add/remove necessary/un-used indices
for ii = 1:numel(coeffs)
	if ~ismember(ii,indexed)
		coeffs{ii}.index = [];
	elseif isempty(coeffs{ii}.index)
		coeffs{ii}.index = (1:size(coeffs{ii}.poly,1)).';
	end
end

% compute the indexed product
integrand = coeffs{1};
for ii = 2:numel(coeffs)
	integrand = (integrand * coeffs{ii});
end
mesh = integrand.mesh;
vals = integrand.vals;
poly = integrand.poly;
indx = integrand.index;

% compute the integral coefficients and collect indexed terms
numDim = size(poly,2) - 1;
intCoeffs = factorial(numDim) * (prod(gamma(1+poly),2) ./ gamma(1+numDim+sum(poly,2)));
intCoeffs = reshape(intCoeffs,1,1,[],1);

[vals,~,indx] = collect(intCoeffs .* vals,[],indx);
vals = vals .* reshape(mesh.volumes,1,1,1,[]);

% distribute indices to the output
if nargout <= 2
	varargout{1} = indx;
else
	varargout = cell(1,min(nargout-1,size(indx,2)));
	for ii = 1:numel(varargout)
		varargout{ii} = indx(:,ii);
	end
end
% DONE
end