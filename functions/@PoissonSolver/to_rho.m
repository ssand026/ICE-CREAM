function [rho] = to_rho(obj,form,data,dim)
% TO_RHO Converts the input data to the weak-formulation charge-density vector
% "rho", where the charge distribution data can be a set of charge-densities,
% wavefunctions, (where the particle's charge is assumed to be +1), or bra-ket
% pairs (used in the calculation of the exact exchange potential). If the data
% contains multiple densities/bras-kets/wavefunctions, the optional argument
% 'dim' can be used to specify the array-dimension/s that separate each set.
arguments (Input)
	obj PoissonSolver
	form (1,1) {mustBeMember(form,["wavefunc","density","bra-ket"])}
	data
	dim (1,:) {mustBeInteger,mustBePositive} = [];
end
arguments (Output)
	rho (:,:) double
end

% validate the inputs
%-----------------------------------------------------------
if (form=="density") || (form=="wavefunc")
	% data must be a double-typed array
	if ~isa(data,"double")
		try data = cast(data,"double");
		catch ME; error("ERROR: the data argument must be a double-typed array")
		end
	end
	if isempty(data)
		error("ERROR: the data array cannot be empty")
	end
	% dim must be a scalar or empty
	if (numel(dim) > 1)
		error("ERROR: the input 'dim' cannot contain more than one element")
	end
elseif (form == "bra-ket")
	% data must be a two-element cell containing double-typed arrays
	if ~isa(data,"cell") || numel(data)~=2
		error("ERROR: the input must be a 2-element cell array")
	end
	if ~isa(data{1},"double")
		try data{1} = cast(data{1},"double");
		catch ME; error("ERROR: the first data value must be a double-typed array")
		end
	end
	if isempty(data{1})
		error("ERROR: the first data value cannot be empty")
	end
	if ~isa(data{2},"double")
		try data{2} = cast(data{2},"double");
		catch ME; error("ERROR: the second data value must be a double-typed array")
		end
	end
	if isempty(data{2})
		error("ERROR: the second data value cannot be empty")
	end
	% dim must have two or fewer elements
	if (numel(dim) > 2)
		error("ERROR: the input 'dim' cannot contain more than two elements")
	end
	% assume sets lie along the same dimension in bra/ket for single input
	if isscalar(dim)
		dim = [dim, dim];
	end
end

% apply the conversion functions
%-----------------------------------------------------------
if (form=="wavefunc")
	% computes rho for a set of wavefunctions
	rho = abs(obj.wavefunc_to_rhs(data,dim));
	
elseif (form=="density")
	% computes rho for a set of charge-densities
	rho = obj.density_to_rhs(data,dim);
	
elseif (form=="bra-ket")
	% computes rho from a single bra-ket pair
	rho = obj.braket_to_rhs(data{1},data{2},dim);
end
% DONE
end