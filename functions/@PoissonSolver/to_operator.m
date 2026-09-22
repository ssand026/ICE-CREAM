function [V] = to_operator(obj,data,dim)
% TO_OPERATOR Converts the input data to the weak-formulation scalar potential
% operator "V". If the data contains multiple potentials, the optional argument
% 'dim' can be used to specify the array-dimension that contains each potential.
arguments
	obj PoissonSolver
	data double
	dim {mustBeInteger,mustBePositive,mustBeScalarOrEmpty} = [];
end

if isempty(dim)
	% single potential specified
	V = obj.scalar_to_operator(data);
else
	% multiple potentials given
	numOperators = size(data,dim);
	V = cell(1,numOperators);
	for ii = 1:numOperators
		V{ii} = obj.scalar_to_operator(part(data,ii,dim));
	end
end
% DONE
end