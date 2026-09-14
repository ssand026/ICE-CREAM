classdef FE_mesh < handle & matlab.mixin.Copyable & matlab.mixin.SetGet
% FE_MESH Handle class for finite-element mesh representation.
%
% FE_MESH is a handle class designed to construct, manage, and process finite-element
% mesh topology, geometric properties, and higher-order polynomial basis nodes across
% arbitrary spatial dimensions. It inherits from handle, matlab.mixin.Copyable, and
% matlab.mixin.SetGet.
%
% The class implements a lazy-evaluation caching system for computationally expensive
% dependent mesh properties—including hypervolumes, barycentric vectors, adjacency matrices,
% boundary faces, nodes, and unit normal vectors. Cached properties are automatically
% invalidated and flushed whenever core mesh coordinates (pts) or connectivity (tri)
% are modified.
%
% FE_MESH also supports higher-order polynomial element meshes. Modifying the polynomial
% degree computes higher-degree nodal placements, sub-simplex midpoint connectivity
% matrices, and polynomial basis exponent arrays required for higher-order finite-element
% basis evaluation.
%
% Syntax
% =======
% obj = FE_MESH(pts,tri)
% obj = FE_MESH(pts,tri,reg)
%
% Input Arguments
% ================
% pts - [numDim x numPts] array describing the mesh node coordinates
% tri - [numVtx x numTri] array describing the mesh connectivity.
% reg - (Optional) [1 x numTri] array used to label different mesh regions.
%
% Name-Value Arguments
% =====================
% degree - Polynomial order for higher-degree finite-element basis functions. Must be a positive
%          integer <= to the number of simplex vertices. Defaults to 1 (linear mesh elements).  
%
% Properties
% ===========
% pts             - [numDim x numPts] array of node coordinates
% tri             - [numVtx x numTri] array containing assigning nodes to each mesh simplex. 
% reg             - Domain/region labels for each simplex.
% degree          - Order of the finite-element representation.
% volumes         - (Dependent) area/volume/hypervolume of each mesh simplex.
% baryvec         - (Dependent) Barycentric vectors associated with each mesh element.
% adjacencyMatrix - (Dependent) Topological node-to-node adjacency matrix.
% boundaryFaces   - (Dependent) List of faces that form the mesh boundary.
% boundaryNodes   - (Dependent) List of node indices that lie on the mesh boundary.
% boundaryNormals - (Dependent) Unit normal vectors computed for the boundary faces.
%
% Pre-Constructor / Static Methods
% =================================
% load - Static pre-constructor method to load and construct a mesh instance from a STL/STEP file or 
%        from one of the many builtin geometry/mesh descriptions.
%
%
% SEE ALSO: FE_COEFF

%% properties
%===========================================================
% core properties
properties (SetObservable)
	pts (:,:) double {mustBeReal,mustBeFinite}
	tri (:,:) uint32 {mustBeInteger,mustBePositive}
	reg (1,:) uint32 {mustBeInteger,mustBeNonnegative}
end
properties (SetObservable,Hidden)
	degree (1,1) uint32 {mustBePositive,mustBeInteger} = 1;
end

% lazily-evaluated properties
properties (Dependent)
	volumes
	baryvec
	adjacencyMatrix
	boundaryFaces
	boundaryNodes
	boundaryNormals
end
% cache stores the values of previously computed dependent properties. This allows changes to
% the mesh to trigger re-computation of dependent properties only when absolutely necessary.
properties (SetAccess=protected,Hidden)
	cache struct = struct;
end

% polynomial-mesh properties
properties (SetAccess=protected,Hidden)
	ptsPoly (:,:) double % polynomial vertex locations
	triPoly (:,:) uint32 % polynomial vertex triangulations
	vtxPoly (:,:) double
	% vtxPoly is a [numPoly x numVtx] array, where vtxPoly(m,n) is nonzero if the polynomial
	% basis-function for the nodes in triPoly(m,:) depend on the barycentric coordinates of the
	% nth simplex vertex. Therefore the polynomial basis-functions for the polynomial vertices
	% in triPoly are given by product  of the barycentric coordinates powers in vtxPoly:
	% basisFunc(triPoly(a,b)) = \prod_n {barycoord(tri(n,b))^vtxPoly(a,n)}
	vtxConn (:,:) logical
	% vtxConn(m,n) is nonzero if the vertex ptsPoly(:,n) sits on the midpoint of the sub-simplex
	% formed by one of the original nodes pts(:,m). For example, if the location of any nodes in
	% the mesh are changed, the location of the polynomial vertex "n" is given by:
	% mean(pts(:,vtxConn(:,n)),2)
end

%% static methods
methods (Static)
	% pre-constructor methods
	%===========================================================
	[varargout] = load(varargin);
	%===========================================================
end

%% constructor method
methods
	%===========================================================
	function obj = FE_mesh(pts,tri,reg,opt)
	% contructs the FE_mesh object
	arguments
		pts (:,:) double {mustBeReal,mustBeFinite,mustBeNonempty}
		tri (:,:) uint32 {mustBeInteger,mustBePositive,mustBeNonempty}
		reg (1,:) uint32 {mustBeInteger,mustBeNonnegative} = [];
		opt.degree (1,1) uint32 {mustBeInteger,mustBePositive} = 1;
	end
	
	% check dim ordering of pts
	if size(pts,1) > size(pts,2)
		pts = pts.';
	end

	% check triangulation limits
	if max(tri,[],"all") > size(pts,2)
		error("ERROR: the given triangulation contains points outside of the node list")
	end

	% check dim ordering of tri
	if (size(tri,1) > size(tri,2))
		if (size(tri,2) == size(pts,1)+1) || (size(tri,2) == size(pts,1)+2)
			tri = tri.';
		else
			error("ERROR: the given triangulation has invalid shape")
		end
	end
	
	% check if region data was appended to the triangulation
	maxVtx = size(pts,1)+1;
	if (size(tri,1) == maxVtx+1)
		warning("the triangulation input to FE_mesh appears to include region data")
		if isempty(reg)
			reg = tri(maxVtx+1,:);
		end
		tri(maxVtx+1:end,:) = [];
	elseif (size(tri,1) > maxVtx+1)
		error("ERROR: the triangulation contains too many vertices per simplex")
	end

	% check region labels
	numTri = size(tri,2);
	if isscalar(reg)
		% do nothing
	elseif isempty(reg)
		% default to region label one 
		reg = 1;
	elseif (numel(reg)~=numTri)
		error("ERROR: the number of region labels should equal the number of triangles")
	end

	% initialize core properties
	obj.pts = pts;
	obj.tri = tri;
	obj.reg = reg;
	
	% initialize degree-dependent properties
	obj.ptsPoly = obj.pts;
	obj.triPoly = obj.tri;
	obj.vtxConn = speye(size(obj.pts,2));
	obj.vtxPoly = eye(size(obj.tri,1));
	
	% initialize polynomial degree
	obj.degree = opt.degree;
	end
	%===========================================================
end

%% misc methods
methods
	%===========================================================
	function [out] = size(obj)
	% returns dimensionality of the mesh, dimensionality of the space,
	% number of nodes, and number of simplexes (mostly for display purposes)
	arguments
		obj FE_mesh
	end
	out = [size(obj.pts,1), size(obj.pts,2), size(obj.tri,1), size(obj.tri,2)];
	end
	%===========================================================
	function [varargout] = num(obj,varargin)
	% Returns the size of the specified mesh parameters
	arguments
		obj FE_mesh
	end
	arguments (Repeating)
		varargin {mustBeMember(varargin,["dim","pts","vtx","tri","poly","reg"])}
	end
	% check mesh degree
	isPolynomial = (obj.degree > 1);
	
	% return the requested mesh-measures
	varargout = cell(1,numel(varargin));
	for ii = 1:numel(varargin)
		switch varargin{ii}
			case "pts"
				if isPolynomial
					varargout{ii} = size(obj.ptsPoly,2);
				else
					varargout{ii} = size(obj.pts,2);
				end
			case "tri"; varargout{ii} = size(obj.tri,2);
			case "dim"; varargout{ii} = size(obj.pts,1);
			case "vtx"; varargout{ii} = size(obj.tri,1);
			case "poly"
				if isPolynomial
					varargout{ii} = size(obj.triPoly,1);
				else
					varargout{ii} = size(obj.tri,1);
				end
			case "reg"; varargout{ii} = max(1,numel(unique(obj.reg)));
			otherwise
				varargout{ii} = NaN;
		end
	end
	% DONE
	end
	%===========================================================
end

%% update triggers for dependent properties
methods
	%===========================================================
	function set.pts(obj,val)
	% clear dependent properties from cache
	dvars = ["volumes","baryvec","boundaryNormals","ptsPoly","triPoly","vtxConn","vtxPoly"];
	for arg = dvars
		if isfield(obj.cache,arg) %#ok
			obj.cache = rmfield(obj.cache,arg); %#ok
		end
	end
	obj.pts = val;
	end
	%===========================================================
	function set.tri(obj,val)
	% clear dependent properties from cache
	dvars = ["volumes","baryvec","adjacencyMatrix",...
		"boundaryFaces","boundaryNodes","boundaryNormals",...
		"ptsPoly","triPoly","vtxConn","vtxPoly"];
	for arg = dvars	
		if isfield(obj.cache,arg) %#ok
			obj.cache = rmfield(obj.cache,arg); %#ok
		end
	end
	obj.tri = val;
	obj.degree = obj.degree; %#ok
	end
	%===========================================================
end

%% evaluate dependent properties
methods
	%===========================================================
	function [out] = get.volumes(obj)
	var = {dbstack().name};
	var = extractAfter(var{1},"get.");
	if ~isfield(obj.cache,var)
		val = hypervolumes(obj.pts,obj.tri);
		obj.cache = setfield(obj.cache,var,val);
	end
	out = getfield(obj.cache,var);
	end
	%===========================================================
	function [out] = get.baryvec(obj)
	var = {dbstack().name};
	var = extractAfter(var{1},"get.");
	if ~isfield(obj.cache,var)
		val = baryvectors(obj.pts,obj.tri);
		val = reshape(val,size(val,1),1,size(val,2),size(val,3));
		obj.cache = setfield(obj.cache,var,val);
	end
	out = getfield(obj.cache,var);
	end
	%===========================================================
	function [out] = get.adjacencyMatrix(obj)
	var = {dbstack().name};
	var = extractAfter(var{1},"get.");
	if ~isfield(obj.cache,var)
		val = adjacency_matrix(obj.tri);
		obj.cache = setfield(obj.cache,var,val);
	end
	out = getfield(obj.cache,var);
	end
	%===========================================================
	function [out] = get.boundaryFaces(obj)
	var = {dbstack().name};
	var = extractAfter(var{1},"get.");
	if ~isfield(obj.cache,var)
		val = boundary_faces(obj.tri);
		obj.cache = setfield(obj.cache,var,val);
	end
	out = getfield(obj.cache,var);
	end
	%===========================================================
	function [out] = get.boundaryNormals(obj)
	var = {dbstack().name};
	var = extractAfter(var{1},"get.");
	if ~isfield(obj.cache,var)
		val = boundary_normals(obj.pts,obj.tri,obj.boundaryFaces);
		obj.cache = setfield(obj.cache,var,val);
	end
	out = getfield(obj.cache,var);
	end
	%===========================================================
	function [out] = get.boundaryNodes(obj)
	var = {dbstack().name};
	var = extractAfter(var{1},"get.");
	if ~isfield(obj.cache,var)
		val = unique(obj.boundaryFaces(:));
		obj.cache = setfield(obj.cache,var,val);
	end
	out = getfield(obj.cache,var);
	end
	%===========================================================
	function [out] = get.ptsPoly(obj)
	var = {dbstack().name};
	var = extractAfter(var{1},"get.");
	if ~isfield(obj.cache,var)
		cache_degree_props(obj);
	end
	out = getfield(obj.cache,var);
	end
	%===========================================================
	function [out] = get.triPoly(obj)
	var = {dbstack().name};
	var = extractAfter(var{1},"get.");
	if ~isfield(obj.cache,var)
		cache_degree_props(obj);
	end
	out = getfield(obj.cache,var);
	end
	%===========================================================
	function [out] = get.vtxConn(obj)
	var = {dbstack().name};
	var = extractAfter(var{1},"get.");
	if ~isfield(obj.cache,var)
		cache_degree_props(obj);
	end
	out = getfield(obj.cache,var);
	end
	%===========================================================
	function [out] = get.vtxPoly(obj)
	var = {dbstack().name};
	var = extractAfter(var{1},"get.");
	if ~isfield(obj.cache,var)
		cache_degree_props(obj);
	end
	out = getfield(obj.cache,var);
	end
	%===========================================================
end

%% supporting methods for polynomial meshes
methods
	%===========================================================
	function set.degree(obj,val)
	% updates the properties associated with polynomial meshes
	
	% check input value
	if (val > obj.num("vtx"))
		error("ERROR: the polynomial order cannot exceed the " + ...
			"number of vertices in the mesh")
	end
	obj.degree = val;
	cache_degree_props(obj)
	end
	%===========================================================
end
methods (Access=private)
	%===========================================================
	function cache_degree_props(obj)
	% updates the polynomial info
	[ptsMid,triMid,ptsInd,triInd] = polymesh_params(obj.pts,obj.tri,obj.degree);
	obj.cache.ptsPoly = ptsMid;
	obj.cache.triPoly = triMid;
	obj.cache.vtxConn = ptsInd;
	obj.cache.vtxPoly = triInd;
	end
	%===========================================================
end
%% DONE
end