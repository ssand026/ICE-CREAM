classdef FE_coeff
% FE_COEFF Data class for scalar, vector, and tensor fields over a mesh
%
% FE_COEFF represents spatially-varying/constant, continuous/disjoint, and scalar/vector/tensor
% fields defined on the supplied finite-element mesh (a FE_mesh object). It allows assembly of the
% weak-form finite-element matrices from the product of arbitrary numers of scalar/vector/tensor
% fields as well as enables one to evaluate gradients, curl, divergence, interpolations, etc.
% of a given field.
%
% Coefficient fields can be initialized via mathematical expressions (function handles or string
% expressions evaluated on the mesh) or constructed directly from existing polynomial basis
% expansion arrays. The class automatically analyzes and tracks field properties including tensor
% rank, maximum polynomial degree, element-wise continuity (continuous vs. disjoint), and spatial
% variation (constant vs. polynomial).
%
% Uniqueness signifiers (index labels) can be attached to individual polynomial terms to govern term
% combination rules during algebraic collection and integration methods. 
%
% This class stores field values internally via their polynomial coefficients. To retrieve the
% real-space node/vertex values, use obj.node_data() or obj.vertex_data() respectively. However, for
% a degree-0 or degree-1 coefficient, the polynomial coefficients are identical to the node values.
%
% Syntax
% =======
% obj = FE_COEFF(mesh,expr)
% obj = FE_COEFF(mesh,expr,poly)
% obj = FE_COEFF(mesh,expr,poly,indx)
%
% Input Arguments
% ================
% mesh - FE_mesh handle object specifying the spatial domain and mesh topology on
%        which the coefficient field is defined.
% expr - Field definition expression or coefficient array. Can be a function handle,
%        string/char evaluation expression, numeric array, or logical mask. When poly
%        is empty or omitted, expr is evaluated on the mesh. When poly is provided,
%        expr must be a numeric array of explicit polynomial coefficient values.
% poly - (Optional) Matrix of size [numTerms x numVtx] containing the powers/exponents
%        of the barycentric coordinates for each term in the polynomial expansion.
%        Defaults to [] when initializing from an expression.
% indx - (Optional) Matrix or vector of positive integer indices specifying uniqueness
%        signifiers for each polynomial coefficient term. Prevents combination of
%        like terms during simplification unless indices match. Defaults to [].
%
% Properties
% ===========
% mesh       - Associated FE_mesh instance on which the coefficient field is defined.
% vals       - Array of numeric coefficients for each term in the polynomial representation.
% poly       - Matrix of powers of barycentric coordinates forming the polynomial representation.
% rank       - Tensor rank of the field (e.g., scalar, vector, or tensor).
% degree     - The overall polynomial degree of the coefficient.
% isDisjoint - Logical flag indicating whether the field is discontinuous across simplex boundaries.
% isConstant - Logical flag indicating whether the field is constant within each individual simplex.
% index      - Uniqueness signifier tags for each polynomial coefficient term.
%
% Pre-Constructor / Static Methods
% =================================
% by_element - pre-constructor method allowing seperate expressions for each vector/tensor element.
% by_region  - re-constructor method allowing seperate expressions for each mesh region.
% integrate  - static method allowing the integration of coefficients over the mesh domain.
%
%
% SEE ALSO: FE_MESH, FE_COEFF.BY_ELEMENT, FE_COEFF.BY_REGION, FE_COEFF.INTEGRATE, FEMAT

%% properties
%===========================================================
properties
	mesh FE_mesh % the mesh on which the coefficient is defined
	vals (:,:,:,:) double % coefficients of each term in the polynomial representation
	poly (:,:) double % powers of the barycoordinates in the polynomial representation
	rank (1,1) double % the rank of the tensor field
	degree (1,1) double % the maximum degree of the polynomial coefficient
	isDisjoint (1,1) logical % whether coefficient is continuous across simplexes
	isConstant (1,1) logical % whether coefficient is constant over each simplex
	index  (:,:) double 
	% index is a uniqueness signifier for each polynomial coefficient, prevents the internal 
	% collection method from combining like polynomials unless they share the same index values,
	% which is required for integration
end

%% static methods
methods (Static)
	%===========================================================
	varargout = by_element(varargin);
	varargout = by_region(varargin);
	varargout = baryfunc(varargin);
	varargout = integrate(varargin);
	%===========================================================
end

%% constructor method
methods
	%===========================================================
	function [obj] = FE_coeff(mesh,expr,poly,indx)
	arguments
		mesh FE_mesh
		expr {mustBeA(expr,{'function_handle','char','string','numeric','logical'})}
		poly (:,:) double {mustBeReal,mustBeFinite} = [];
		indx (:,:) double {mustBeInteger,mustBePositive} = [];
	end
	
	if (nargin <= 2)
		% initialize a FE_coeff, convert the expression to polynomial coefficients
		%-----------------------------------------------------------
		% evaluate the expression
		data = mesh.evalfun(expr);
		% check the size/ordering of the mesh data
		[data,isDisjoint,isConstant,rank] = check(mesh,data);
		% convert data to vertex-data
		if isConstant
			% data is domain/simplex constant
			poly = zeros(1,size(mesh.vtxPoly,2));
		else
			poly = mesh.vtxPoly;
			if ~isDisjoint
				% data is defined over nodes
				data = nodal_to_vertex(mesh,data);
			end
		end
		% convert to polynomial basis
		vals = vertex_to_poly(data,poly);
	else
		% reconstruct a FE_coeff using existing polynomial coefficients
		%-----------------------------------------------------------
		if ~isnumeric(expr)
			error("ERROR: the second input must be a numeric array corresponding" + ...
				"to the polynomial coefficients when the polynomial powers are given.")
		else
			vals = double(expr);
		end
		% ensure index list stores simplex labels along each column
		if isvector(indx)
			indx = indx(:);
		end
		% check the size/ordering of data
		[vals,isDisjoint,isConstant,rank] = check(mesh,vals);
	end
	
	% assemble the FE_coeff object
	obj.mesh = mesh;
	obj.vals = vals;
	obj.poly = poly;
	obj.index = indx;
	obj.rank = rank;
	obj.degree = max(sum(poly,2));
	obj.isDisjoint = isDisjoint;
	obj.isConstant = isConstant;
	% DONE
	end
	%===========================================================
end

%% operator overload
methods
	%===========================================================
	function [out] = size(obj,dims)
	arguments
		obj FE_coeff
		dims (1,:) double {mustBeInteger,mustBePositive} = [1:ndims(obj.vals)];
	end
	out = size(obj.vals,dims);
	end
	%===========================================================
	function obj = transpose(obj)
	obj.vals = permute(obj.vals,[2,1,3:ndims(obj.vals)]);
	end
	%===========================================================
	function obj = ctranspose(obj)
	obj.vals = permute(obj.vals,[2,1,3:ndims(obj.vals)]);
	if ~isreal(obj.vals)
		obj.vals = conj(obj.vals);
	end
	end
	%===========================================================
	function [out] = double(obj)
	out = obj.vals;
	end
	%===========================================================
end

%% other methods
methods
	%===========================================================
	function [out] = vertex_data(obj)
	% return the coefficient's value at every mesh point
	if obj.isConstant && ~obj.isDisjoint
		% coefficient is constant over the domain
		out = repmat(obj.vals,1,1,obj.mesh.num("poly"),obj.mesh.num("tri"));
		return
	elseif obj.isConstant && obj.isDisjoint
		% coefficient is constant over each simplex
		out = poly_to_vertex(obj.vals,obj.poly);
		out = repmat(out,1,1,obj.mesh.num("poly"));
	else
		% coefficient is non-constant
		[out] = poly_to_vertex(obj.vals,obj.poly);
	end
	end
	%===========================================================
	function [out] = node_data(obj)
	% return the coefficient's value at every mesh node
	if obj.isConstant && ~obj.isDisjoint
		% coefficient is constant over the domain
		out = repmat(obj.vals,1,1,obj.mesh.num("pts"));
		return
	elseif obj.isConstant && obj.isDisjoint
		% coefficient is constant over each simplex
		vtxData = poly_to_vertex(obj.vals,obj.poly);
		vtxData = repmat(vtxData,1,1,obj.mesh.num("poly"));
		[out] = vertex_to_nodal(obj.mesh,vtxData);
	else
		% coefficient is non-constant
		[vtxData] = poly_to_vertex(obj.vals,obj.poly);
		[out] = vertex_to_nodal(obj.mesh,vtxData);
	end
	end
	%===========================================================
	function [obj] = collected(obj)
	% combine coefficient terms with the same polynomial components
	[obj.vals,obj.poly] = collect(obj.vals,obj.poly);
	[~,obj.isDisjoint,obj.isConstant] = check(obj.mesh,obj.vals);
	end
	%===========================================================
	function [obj] = expanded(obj)
	% expands a coefficient so its value is defined on each simplex
	obj.vals = expand(obj.mesh,obj.vals);
	[~,obj.isDisjoint,obj.isConstant] = check(obj.mesh,obj.vals);
	end
	%===========================================================
end
%% DONE
end