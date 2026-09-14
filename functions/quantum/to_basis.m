function [data] = to_basis(data,target,intoB,exitB)
% TO_BASIS: transfers data in/out of the subspace of some larger/complete basis
% 
% SYNTAX: 
% [out] = TO_BASIS(data,target,intoB,exitB)
%
% The input 'data' can take the form of an operator or a set of row/column-vectors.
% The form of the data is determined automatically from the full/reduced basis sizes,
% but in general, square matrices will be interpreted as operators, and
% non-square inputs are interpreted as row/column-vectors.
% 
% The input 'target' is used to denote the target basis for the transformed data.
% Specifying "full" will return the complete-basis representation, while setting
% 'target' equal to "reduced" will return the subspace form of the data.
% 
% For a basis of length 'n', and subspace of length 'm':
%   'intoB' is a [m x n] reduction matrix that maps the full basis to its subspace
%   'exitB' is a [n x m] expansion matrix that maps the subspace to the full basis
% These matrices should satsify the relation: (intoB * exitB) = eye(m,m). While only
% the matrix 'intoB' OR 'exitB" is required, omitting the conversion matrix needed for
% the transformation will result in the use of more expensive linear solves.
%
% SEE ALSO: FEMAT_DIRICHLET, FEMAT_NEUMANN
arguments
	data {mustBeA(data,{'numeric','cell'})}
	target (1,1) {mustBeMember(target,["full","reduced"])}
	intoB (:,:) double = [];
	exitB (:,:) double = [];
end

%% check inputs
%===========================================================
% check for conversion matrix
if isempty(intoB) && isempty(exitB)
	error("ERROR: neither the regular nor inverted basis conversion matrix was given.")
end

% check intoB
if ~isempty(intoB)
	% ensure rows(intoB) < cols(intoB)
	if size(intoB,1) >= size(intoB,2)
		error("ERROR: the basis conversion matrix must have more columns than rows.")
	end
	
	if isempty(exitB)
		% check if intoB is a pseudo-permutation matrix
		nonZeros = nonzeros(intoB);
		if isempty(nonZeros)
			error("ERROR: basis matrix is all zero-valued.");
		end
		a = nonZeros(1);
		if all(nonZeros==a) && all(full(sum(intoB,2))==a)
			% exitB = inv(a*intoB) = transpose(a*intoB)/(a*a)
			exitB = intoB.'/(a^2);
		end
	end
end

% check exitB
if ~isempty(exitB)
	% ensure rows(exitB) > cols(exitB)
	if (size(exitB,1) <= size(exitB,2))
		error("ERROR: the inverted basis conversion matrix must have more rows than columns.")
	end
	
	if isempty(intoB)
		% check if exitB is a pseudo-permutation matrix
		nonZeros = nonzeros(exitB);
		if isempty(nonZeros)
			error("ERROR: inverted basis matrix is all zero-valued.");
		end
		a = nonZeros(1);
		if all(nonZeros==a) && all(full(sum(exitB,1))==a)
			% intoB = inv(a*exitB) = transpose(a*exitB)/(a*a)
			intoB = exitB.'/(a^2);
		end
	end
end

% check that the sizes of intoB and exitB agree
if ~isempty(intoB) && ~isempty(exitB)
	if ~isequal(size(intoB,[1,2]),size(exitB,[2,1]))
		error("ERROR: the sizes of the basis matrices are incompatible.")
	end
	[minSize,maxSize] = size(intoB);
elseif ~isempty(intoB)
	[minSize,maxSize] = size(intoB);
elseif ~isempty(exitB)
	[maxSize,minSize] = size(exitB);
end

%% convert data
%===========================================================
if iscell(data)
	isvalid = cellfun(@isnumeric,data,"UniformOutput",true);
	if ~all(isvalid)
		error("ERROR: the contents of the cell array must be numeric.")
	end
	for ii = 1:numel(data)
		if isempty(data{ii}); continue; end
		data{ii} = convert(double(data{ii}));
	end
else
	data = convert(double(data));
end

%% internal functions
%***********************************************************
function [datum] = convert(datum)
%% check the size/shape of the data
%===========================================================
switch size(datum,1)
	case minSize; rowSize = 1;
	case maxSize; rowSize = 2;
	otherwise; rowSize = 0;
end

switch size(datum,2)
	case minSize; colSize = 1;
	case maxSize; colSize = 2;
	otherwise; colSize = 0;
end

if (rowSize==1) && (colSize==1)
	% reduced basis matrix
	isFull = false;
	shape = "mat";
elseif (rowSize==1) && (colSize==0)
	% reduced basis row-vector
	isFull = false;
	shape = "row";
elseif (rowSize==0) && (colSize==1)
	% reduced basis column-vector
	isFull = false;
	shape = "col";
elseif (rowSize==2) && (colSize==2)
	% full basis matrix
	isFull = true;
	shape = "mat";
elseif (rowSize==2) && (colSize==0)
	% full basis row-vector
	isFull = true;
	shape = "row";
elseif (rowSize==0) && (colSize==2)
	% full basis column-vector
	isFull = true;
	shape = "col";
elseif (rowSize==0) && (colSize==0)
	% no overlap with the transformation matrix
	error("ERROR: the size of the data is incompatible with the size of the conversion matrix.")
elseif (target == "full")
	% assume data is in the reduced basis
	isFull = false;
	if (rowSize==2) && (colSize==1);  shape = "col"; end
	if (rowSize==1) && (colSize==2);  shape = "row"; end
	
elseif (target == "reduced")
	% assume data is in the full basis
	isFull = true;
	if (rowSize==2) && (colSize==1); shape = "row"; end
	if (rowSize==1) && (colSize==2); shape = "col"; end
end

% check if data is already in the correct basis
if (isFull && (target=="full")) || (~isFull && (target=="reduced"))
	warning("The data is already in the "+target+" basis.")
	return
end

%% convert data to the proper basis
%===========================================================
if (shape == "mat")
	% data is an operator
	%-----------------------------------------------------------
	% check if data is a matrix
	if ~ismatrix(datum)
		error("ERROR: the input data cannot be multi-dimensional " + ...
			"when the first two dimensions have equal length.")
	end
	
	% check the matrix symmetry
	isHerm = ishermitian(datum) - ishermitian(datum,"skew");
	isSymm = issymmetric(datum) - issymmetric(datum,"skew");
	
	if isFull
		% convert to reduced basis
		if ~isempty(intoB)
			datum = (intoB * datum * intoB.');
		else
			datum = (exitB \ datum / exitB.');
		end
	else
		% convert to full basis
		if ~isempty(exitB)
			datum = (exitB * datum * exitB.');
		else
			datum = (intoB \ datum / intoB.');
		end
	end
	
	% restore the matrix symmetry
	if (isHerm == 0) && (isSymm == 0)
		% no symmetries
	elseif (isHerm == -1)
		% skew hermitian
		datum = (datum - datum')/2;
	elseif (isHerm == +1)
		% hermitian
		datum = (datum + datum')/2;
	elseif (isSymm == -1)
		% skew symmetric
		datum = (datum - datum.')/2;
	elseif (isSymm == +1)
		% symmetric
		datum = (datum + datum.')/2;
	end
else
	% data is a set of row/column vectors
	%-----------------------------------------------------------
	% move basis-commuting dimension to the front
	if (shape == "row"); dimPerm = [1,2,3:ndims(datum)]; end
	if (shape == "col"); dimPerm = [2,1,3:ndims(datum)]; end
	datum = permute(datum,dimPerm);
	otherDims = size(datum,[2:numel(dimPerm)]);
	
	if isFull
		% convert to reduced basis
		datum = reshape(datum,maxSize,[]);
		if ~isempty(intoB)
			datum = intoB * datum;
		else
			datum = exitB \ datum;
		end
		datum = reshape(datum,[minSize,otherDims]);
	else
		% convert to full basis
		datum = reshape(datum,minSize,[]);
		if ~isempty(exitB)
			datum = exitB * datum;
		else
			datum = intoB \ datum;
		end
		datum = reshape(datum,[maxSize,otherDims]);
	end
	
	% restore the original dimension order
	datum = permute(datum,dimPerm);
end

end
%***********************************************************
%% DONE
end