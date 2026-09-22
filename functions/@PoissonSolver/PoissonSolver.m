classdef PoissonSolver < handle & matlab.mixin.Copyable & matlab.mixin.SetGet
% POISSONSOLVER is a handle class that constructs and solves finite-element representations of the
% inhomogeneous/screened Poisson equation:
%
%      -∇·(ε(r) ∇ V(r)) + κ²(r)V(r) = ρ(r)
% where:
%   ε is the permittivity
%   κ² is the screening coefficient
%   ρ is the charge density
%   V is the electrostatic potential
%
% Boundary conditions are established using a truncated multipole expansion. By embedding
% the primary interior mesh within a spherical boundary domain where charge density is
% zero and material properties are constant, higher-order multipole expansion terms
% integrate to zero on the outer boundary surface, yielding the constraint:
%      V_tot = Q / (4*pi*ε*R)
% where:
%   V_tot is the integral of the potential V over the spherical boundary
%   Q is the total charge enclosed within the sphere
%   ε is the permittivity on the spherical boundary
%   R is the radius of the enclosing sphere
%
% The class supports inhomogeneous and anisotropic permittivity, inhomogeneous screening
% coefficients, and preconditioning of the linear system/factorizations for fast repeated solves.
%
% REQUIRES the matlab package GIBBON to generate the enclosing sphere meshes for 3D systems.
%
% Syntax
% =======
% obj = POISSONSOLVER(intMesh)
% obj = POISSONSOLVER(intMesh,intPermit,intScreen)
% obj = POISSONSOLVER(intMesh,intPermit,intScreen,extPermit,extScreen)
% obj = POISSONSOLVER(___,Name,Value)
%
% Input Arguments
% ================
% intMesh      - FE_mesh object defining the interior primary domain.
% intPermit    - (Optional) Interior domain permittivity coefficient. Can be an FE_coeff
%                object, function handle, char/string expression, numeric array, or logical
%                mask. Defaults to vacuum permittivity (ε0) if empty or omitted.
% intScreen    - (Optional) Interior domain Debye-Hückel screening coefficient (κ²).
%                Can be an FE_coeff object, function handle, string expression, numeric
%                array, or logical mask. Defaults to 0 if empty or omitted.
% extPermit    - (Optional) Scalar numeric permittivity constant for the exterior domain
%                outside the primary mesh. Defaults to vacuum permittivity (ε0) if empty/omitted.
% extScreen    - (Optional) Scalar numeric screening constant for the exterior domain
%                outside the primary mesh. Defaults to 0 if empty or omitted.
%
% Name-Value Arguments
% =====================
% method       - Solution algorithm / preconditioner method for linear system solves. Options
%                include "decompose", "invert", "lsqr", "gmres", "bicg", "bicgstab", "bicgstabl",
%                "cgs", "qmr", "tfqmr", "minres", "symmlq", or "pcg". Defaults to "decompose".
% ordering     - Matrix reordering strategy to minimize fill-in during matrix factorization.
%                Options are "amd", "symamd", "symrcm", or "dissect". Defaults to "dissect".
% exteriorMesh - User-supplied FE_mesh object specifying a pre-built spherical exterior domain.
% radiusScale  - Scale factor for the enclosing sphere radius relative to the interior domain.
%                Must be >= 1. Defaults to 1.05.
% exteriorGap  - Gap factor controlling the exterior mesh domain thickness. Must be >= 1.
%                Defaults to 1.5.
%
% Properties
% ===========
% intMesh     - FE_mesh handle object defining the interior domain.
% extMesh     - FE_mesh handle object defining the combined exterior spherical domain.
% intSize     - Scalar numeric value specifying the number of nodes in the interior mesh.
% extSize     - Scalar numeric value specifying the total number of nodes in the exterior mesh.
% L           - Reordered finite-element assembly matrix for the Poisson system.
% dL          - Preconditioned matrix solver structure enabling accelerated linear solves.
% isScreened  - Logical flag indicating whether Debye-Hückel screening is active.
% solveOrder  - Vector permutation indices applied to reorder the system matrix.
% basisOrder  - Inverse permutation indices used to map system solution back to original node ordering.
% method      - Active linear system solution / preconditioner method name.
%
% Pre-Constructor / Static Methods
% =================================
% init_to_rhs      - Static helper method initializing the right-hand side conversion functions.
% init_to_operator - Static helper method initializing potential operator generator functions.
%
%
% SEE ALSO: FE_MESH, FE_COEFF, ENCLOSING_SPHERE, PRECONDITIONED

%% properties
%===========================================================
properties (SetAccess=protected)
	intMesh FE_mesh
	extMesh FE_mesh
	intSize (1,1) double
	extSize (1,1) double
	L double {mustBeMatrix}
	dL preconditioned
	isScreened (1,1) logical
	solveOrder (1,:) double
	basisOrder (1,:) double
end

properties
	method (1,1) {mustBeMember(method, ["decompose","invert","iter", ...
		"lsqr","gmres","bicg","bicgstab","bicgstabl","cgs","qmr","tfqmr", ...
		"minres","symmlq","pcg"])} = "decompose";
end

% pre-processed conversion functions
properties (Hidden,SetAccess=protected)
	wavefunc_to_rhs function_handle
	density_to_rhs function_handle
	braket_to_rhs function_handle
	scalar_to_operator function_handle
end

%% static methods
methods (Static)
	% initializers for the conversion functions
	%===========================================================
	[varargout] = init_to_rhs(varargin);
	[varargout] = init_to_operator(varargin);
	%===========================================================
end

%% constructor
methods
	%===========================================================
	function [obj] = PoissonSolver(intMesh,intPermit,intScreen,extPermit,extScreen,opt)
	% Returns a function that gives solutions to the inhomogeneous/screened Poisson equation
	%
	% Input Arguments
	% ================
	%	intMesh: an FE_mesh object
	%   intPermit: the permittivity coefficient in the mesh interior
	%   intScreen: the screening coefficient in the mesh interior
	%   extPermit: the permitivity constant in the region outside the mesh
	%   extScreen: the screening constant in the region outside the mesh
	%
	% The screening coefficients refer to the squared inverse of the debeye length
	%
	% SEE ALSO: ENCLOSING_SPHERE, INIT_TO_RHS, INIT_TO_OPERATOR
	arguments
		intMesh   FE_mesh
		intPermit {mustBeA(intPermit,{'FE_coeff','function_handle','char','string','numeric','logical'})} = [];
		intScreen {mustBeA(intScreen,{'FE_coeff','function_handle','char','string','numeric','logical'})} = [];
		extPermit double {mustBeScalarOrEmpty,mustBeFinite} = [];
		extScreen double {mustBeScalarOrEmpty,mustBeFinite} = [];
		% optional arguments
		opt.method (1,1) {mustBeMember(opt.method,["decompose","invert","lsqr","gmres", ...
			"bicg","bicgstab","bicgstabl","cgs","qmr","tfqmr","minres","symmlq","pcg"])} = "decompose";
		opt.ordering (1,1) {mustBeMember(opt.ordering,["amd","symamd","symrcm","dissect"])} = "dissect";
		% exterior mesh options
		opt.exteriorMesh FE_mesh
		opt.radiusScale (1,1) double {mustBeGreaterThanOrEqual(opt.radiusScale,1)} = 1.05;
		opt.exteriorGap (1,1) double {mustBeGreaterThanOrEqual(opt.exteriorGap,1)} = 1.5;
	end
	
	% set defualt values for the permittivity/screening coefficients
	if isempty(intPermit); intPermit = atomic_constants("ε0"); end
	if isempty(extPermit); extPermit = atomic_constants("ε0"); end
	if isempty(intScreen); intScreen = 0; end
	if isempty(extScreen); extScreen = 0; end
	
	
	% attempt evaluation of interior coefficients and extract vertex data
	%-----------------------------------------------------------
	if ~isa(intPermit,"FE_coeff")
		intPermit = FE_coeff(intMesh,intPermit);
	end
	if ~isa(intScreen,"FE_coeff")
		intScreen = FE_coeff(intMesh,intScreen);
	end
	intPermit = vertex_data(intPermit);
	intScreen = vertex_data(intScreen);
	
	% prepare the spherical exterior mesh
	%-----------------------------------------------------------
	if isfield(opt,"exteriorMesh")
		% load in the enclosing mesh
		extMesh = opt.exteriorMesh;
		
		% get the mesh radius
		bndPts = extMesh.pts(:,extMesh.boundaryNodes);
		radii = vecnorm(bndPts-mean(bndPts,2),2,1);
		R = mean(radii);
		
		% check the mesh validity
		if (std(radii) > 0.01*R)
			error("ERROR: the exterior mesh appears to have a non-spherical boundary")
		end
	else
		% Construct a mesh where the original domain is enclosed in a solid sphere
		opts = {"exteriorGap",opt.exteriorGap,"radiusScale",opt.radiusScale};
		[~,R,extPts,extTri,extReg] = enclosing_sphere(intMesh.pts,intMesh.tri,intMesh.reg,opts{:});
		extMesh = FE_mesh(extPts,extTri,extReg);
	end
	numIntPts = intMesh.num("pts");
	numExtPts = extMesh.num("pts");
	
	% check if screening effects are enabled
	%-----------------------------------------------------------
	if isequal(extScreen,0) && all(intScreen(:)==0)
		isScreened = false;
	else
		isScreened = true;
	end
	
	% assemble the combined Poisson matrix
	%-----------------------------------------------------------
	% get the region labels for the interior/exterior regions
	extLabel = max(extMesh.reg);
	intLabel = setdiff(unique(extMesh.reg),extLabel);
	
	% incorporate permittivity
	permCoeff = FE_coeff.by_region(extMesh,intLabel,intPermit,extLabel,extPermit);
	L = FEmat(extMesh,"stiffness",permCoeff);
	
	% incorporate screening
	if isScreened
		screenCoeff = FE_coeff.by_region(extMesh,intLabel,intScreen,extLabel,extScreen);
		L = L + FEmat(extMesh,"scalar",screenCoeff);
	end
	
	% apply Lagrange multipliers to enforce the multipole-expansion constraint
	%-----------------------------------------------------------
	% If the relative permitivity 'perm' is constant on a spherical surface of
	% radius 'R' enclosing a charge 'Q', the average potential on the surface
	% of the sphere must equal Q/(4*pi*perm*R),
	
	% compute the surface integral vector (summed fractional surface area of each node node)
	faces = extMesh.boundaryFaces;
	faceAreas = hypervolumes(extMesh.pts,faces.');
	faceAreas = repmat(faceAreas(:),1,size(faces,2));
	surfIntegral = sparse(1,double(faces),faceAreas,1,numExtPts);
	surfIntegral = surfIntegral / sum(surfIntegral);
	
	% apply the boundary-condition/Lagrange-multiplier to L
	bc = (4*pi*extPermit*R) * surfIntegral; % (4*pi*ε*R) * surface-integral vector
	if isScreened
		% settle for a QR solver until I can verify that the dual Lagrange-multipliers
		% work when L is not a stiffness matrix
		L(end+1,:) = bc;
	else
		% Add a extra degree of freedom to the potential so that the Poisson
		% matrix "L" is square and symmetric, as opposed to rectangular. This
		% enables the use of faster solve routines i.e. (LDL decomp vs. QR decomp).
		% However, v = L\rhs will now contain extra rows, which must be discarded.
		a = mean(nonzeros(bc));
		L =[[L, bc', bc'];
			[bc, -a, +a];
			[bc, +a, -a]];
	end
	
	% re-order the matrix to minimize the factorization-size of L
	%-----------------------------------------------------------
	% apply a strict amd reordering
	spvals = spparms;
	spparms("tight");
	switch opt.ordering
		case "amd"
			% use the Approximate Minimum Degree ordering (best for iterative solvers)
			obj.solveOrder = amd(L);
		case "symamd"
			% use the symmetric Approximate Mininimum Degree ordering
			obj.solveOrder = symamd(L);
		case "symrcm"
			% use the symmetric Reverse Cuthill-McKee ordering algorithm
			obj.solveOrder = symrcm(L);
		case "dissect"
			% use the dissection ordering algorithm (best for the default '\' solver)
			obj.solveOrder = dissect(L,"NumIterations",32);
	end
	spparms(spvals);
	
	% re-sort the rows/columns of L
	[~,obj.basisOrder] = sort(obj.solveOrder);
	L = L(obj.solveOrder,obj.solveOrder);
	
	% compute any pre-conditioners for the matrix to speed up repeated linear solves
	dL = preconditioned(L,"method",opt.method,"warnings",true);
	dL.warnings = false;
	
	% populate the class properties
	%-----------------------------------------------------------
	obj.L = L;
	obj.dL = dL;
	obj.method = dL.method; % use the validated preconditioner method
	obj.intMesh = intMesh;
	obj.extMesh = extMesh;
	obj.intSize = numIntPts;
	obj.extSize = numExtPts;
	obj.isScreened = isScreened;
	
	% initialize the conversion functions
	%-----------------------------------------------------------
	[func1,func2,func3] = PoissonSolver.init_to_rhs(intMesh);
	[func4] = PoissonSolver.init_to_operator(intMesh);
	obj.wavefunc_to_rhs = func1;
	obj.density_to_rhs = func2;
	obj.braket_to_rhs = func3;
	obj.scalar_to_operator = func4;
	end
	%===========================================================
end

%% core methods
methods
	%===========================================================
	function [v_int,v_ext] = solve(obj,rhs)
	% SOLVE Solves for the Poisson potential given the right-side vector "rhs"
	arguments (Input)
		obj PoissonSolver
		rhs (:,:) double
	end
	% get the size of the internal/external mesh
	numIntPts = obj.intSize;
	numExtPts = obj.extSize;
	
	% set charge to zero outside the interior region
	rhs([(numIntPts+1):numExtPts],:) = 0;
	
	% append the total enclosed charge "Q" to rhs to enforce the boundary condition
	if obj.isScreened
		rhs(end+1,:) = sum(rhs,1);
	else
		rhs(end+(1:2),:) = sum(rhs,1);
	end
	
	% update the solver method and return the potential
	rhs = rhs(obj.solveOrder,:);
	obj.dL.method = obj.method;
	obj.method = obj.dL.method;
	v_ext = (obj.dL \ rhs);
	v_ext = v_ext(obj.basisOrder,:);
	
	% return the potentials for the interior/exterior regions
	v_ext = v_ext(1:numExtPts,:);
	v_int = v_ext(1:numIntPts,:);
	end
	%===========================================================
	function [rho] = isolve(obj,v)
	% ISOLVE Solves for the charge density "rho" given a potential "v"
	arguments (Input)
		obj PoissonSolver
		v (:,:) double
	end
	
	numIntPts = obj.intSize;
	numExtPts = obj.extSize;
	LL = obj.L(obj.basisOrder,obj.basisOrder);
	
	if (size(v,1) ~= numIntPts) && (size(v,1) ~= numExtPts)
		error("ERROR: the number of rows in the potential must be compatible with the mesh size")
	elseif size(v,1) == numIntPts
		rho = LL(1:numIntPts,1:numIntPts) * v;
	elseif size(v,1) == numExtPts
		rho = LL(1:numExtPts,1:numExtPts) * v;
	end
	end
	%===========================================================
end
end