function [obj] = by_region(mesh,regIndex,expr,opt)
% constructs a FE_coeff from a list of regions and the coefficient values/expressions corresponding
% to each region 
arguments
	mesh FE_mesh
end
arguments (Repeating)
	regIndex (1,:) {mustBeInteger,mustBeNonnegative}
	expr {mustBeA(expr,{'function_handle','char','string','numeric','logical'})}
end
arguments
	opt.defaultTo {mustBeA(opt.defaultTo,{'function_handle','char','string','numeric'})} = [];
end

% fill unspecified region values using the given default
if ~isempty(opt.defaultTo)
	in_args = unique([regIndex{:}]);
	in_mesh = unique(mesh.reg(:));
	unspec = setdiff(in_mesh,in_args);
	if ~isempty(unspec)
		regIndex{end+1} = unspec(:).';
		expr{end+1} = opt.defaultTo;
	end
end

% check region indices
in_args = unique([regIndex{:}]);
in_mesh = unique(mesh.reg(:));
unspec = setdiff(in_mesh,in_args); % regions without defined values
exspec = setdiff(in_args,in_mesh); % regions not in mesh
if ~isempty(unspec)
	invalid = join(""+unspec,",");
	error("ERROR: values not defined for mesh regions ["+invalid+"]");
end
if ~isempty(exspec)
	invalid = "["+join(""+exspec,",")+"]";
	error("ERROR: listed regions "+invalid+" are not found within the mesh");
end

% compute the expression values in each region
for ii = 1:numel(expr)
	mesh_ii = mesh.subset("reg",regIndex{ii});
	expr{ii} = mesh_ii.evalfun(expr{ii});
	expr{ii} = expand(mesh_ii,expr{ii});
end

% ensure that the tensor-rank in every region is compatible
[numDim] = mesh.num("dim");
isScalar = cellfun(@(x)isequal(size(x,[1,2]),[1,1]),expr);
isMatrix = cellfun(@(x)isequal(size(x,[1,2]),[numDim,numDim]),expr);
isVector = ~(isScalar | isMatrix);
if any(isVector) && (any(isScalar) || any(isMatrix))
	error("ERROR: the rank of the coefficients are not the same across regions ")
elseif (any(isScalar) && any(isMatrix))
	for ii = 1:numel(expr)
		% promote to multiple of the identity matrix
		if isScalar(ii)
			expr{ii} = eye(numDim) .* expr{ii};
		end
	end
end

% combine region data into a single coefficient
data = zeros(size(expr{1}));
for ii = 1:numel(expr)
	is_reg_ii = ismember(mesh.reg,regIndex{ii});
	data(:,:,:,is_reg_ii) = expr{ii};
end

% invoke constructor
[obj] = FE_coeff(mesh,data);

% DONE
end