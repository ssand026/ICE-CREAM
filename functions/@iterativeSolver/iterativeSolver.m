classdef iterativeSolver < handle & matlab.mixin.Copyable & matlab.mixin.SetGet
% ITERATIVESOLVER Iterative linear solver with automatic method selection.
%
% OBJ = ITERATIVESOLVER(A) constructs a solver object for repeatedly
% solving linear systems involving the matrix A.
% 
% OBJ = ITERATIVESOLVER(A,METHOD) requests a specific iterative solver.
% If the requested method is incompatible with the matrix properties,
% a suitable alternative is selected automatically.
% 
% The solver supports both square and rectangular matrices and provides
% optional matrix equilibration and incomplete-factorization
% preconditioners to accelerate convergence.
% 
% Supported iterative methods
% ===========================
%   lsqr       Least-squares QR iteration
%   gmres      Generalized minimum residual *
%   bicg       Bi-conjugate gradient
%   bicgstab   Stabilized bi-conjugate gradient
%   bicgstabl  BiCGStab(l)
%   cgs        Conjugate gradient squared *
%   qmr        Quasi-minimal residual
%   tfqmr      Transpose-free QMR
%   minres     Minimum residual
%   symmlq     Symmetric LQ
%   pcg        Preconditioned conjugate gradient *
% 
% Matrix analysis
% ---------------
% During construction the matrix is analyzed to determine whether it is:
%   - Sparse or dense
%   - Square or rectangular
%   - Symmetric
%   - Hermitian
%   - Positive definite
% 
% These properties are used to validate solver choices and configure
% compatible preconditioners.
% 
% Configurable Properties
% =======================
% equilibrate:
%	If set to true, apply row/column scaling using EQUILIBRATE.
% factorize:
%	Construct incomplete LU or incomplete Cholesky preconditioners when applicable.
% useCholFact 
%	Use ICHOL for positive-definite systems solved with PCG, MINRES, or SYMMLQ.
% checkPosDef:
%	Explicitly test positive definiteness using a Cholesky factorization. If 
%	disabled, symmetric and Hermitian matrices are assumed to be positive definite.
% parallel:
%	Solve multiple right-hand sides in parallel when supported by the selected method.
%
% SEE ALSO:
%	ITERATIVESOLVER.SOLVE, ITERATIVESOLVER.LITSOLVE
%	EQUILIBRATE, ILU, ICHOL,
%	LSQR, GMRES, BICG, BICGSTAB, BICGSTABL, CGS, QMR, TFQMR, MINRES, SYMMLQ, PCG


%% properties
properties (SetAccess=protected)
	A (:,:) double {mustBeMatrix}
	isSparse (1,1) logical = false;
	isSquare (1,1) logical = false;
	isSymm (1,1) logical = false;
	isHerm (1,1) logical = false;
	isPosDef (1,1) logical = false;
end

properties (SetObservable)
	method (1,1) string = "auto";
end

% solver settings/options
properties
	tolerance (1,1) double {mustBeInRange(tolerance,0,1)} = 1e-12;
	maxIters (1,1) uint16 {mustBePositive,mustBeInteger} = 256;
	restarts (1,1) uint16 {mustBeInteger} = 0;
	equilibrate (1,1) logical = false;
	factorize (1,1) logical = true;
	useCholFact (1,1) logical = false;
	checkPosDef (1,1) logical = false;
	warnings (1,1) logical = true;
	parallel (1,1) logical = false;
end

% preconditioners
properties (Hidden)
	R (:,1) double {mustBeVector} = []; % row-scaling vector from equilibrate
	C (:,1) double {mustBeVector} = []; % column-scaling vector from equilibrate
	B (:,:) double {mustBeMatrix} = []; % equilibrated matrix (R .* A .* C')
	L (:,:) double {mustBeMatrix} = []; % lower ilu component
	U (:,:) double {mustBeMatrix} = []; % upper ilu component
	Q (:,:) double {mustBeMatrix} = []; % lower ichol component
end

%% Static methods
methods (Static)
	%===========================================================
	[varargout] = check_method(varargin);
	[varargout] = litsolve(varargin);
	%===========================================================
end

%% constructor method
methods
	%===========================================================
	function [obj] = iterativeSolver(A,arg,opt)
	arguments
		A (:,:) double {mustBeMatrix}
		arg.method (1,1) string 
		opt.tolerance (1,1) double {mustBeInRange(opt.tolerance,0,1)}
		opt.maxIters (1,1) uint16 {mustBePositive}
		opt.restarts (1,1) uint16
		opt.equilibrate (1,1) logical
		opt.factorize   (1,1) logical
		opt.useCholFact (1,1) logical
		opt.checkPosDef (1,1) logical
		opt.warnings (1,1) logical
		opt.parallel (1,1) logical
	end
	
	% overwrite any pre-set optional properties 
	set(obj,opt);

	% evaluate matrix properties
	%-----------------------------------------------------------	
	if obj.checkPosDef
		[isSquare,isSymm,isHerm,isPosDef] = matrix_properties(A);
	else
		% assume positive-definiteness if matrix is symmetric/hermitian
		[isSquare,isSymm,isHerm] = matrix_properties(A);
		isPosDef = (isSymm || isHerm);
	end

	obj.A = A;
	obj.isSparse = issparse(A);
	obj.isSquare = isSquare;
	obj.isSymm = isSymm;
	obj.isHerm = isHerm;
	obj.isPosDef = isPosDef;

	
	% disable invalid conditioners
	%-----------------------------------------------------------
	if (obj.equilibrate == true) && ~obj.isSquare
		obj.equilibrate = false;
		if obj.warnings
			warning("Disabling equilibrate preconditioning, matrix is non-square.");
		end
	end
	if (obj.equilibrate == true) && (obj.isSymm || obj.isHerm)
		obj.equilibrate = false;
		if obj.warnings
			warning("Disabling equilibrate preconditioning, matrix is symmetric.");
		end
	end
	if (obj.factorize == true) && ~obj.isSparse
		obj.factorize = false;
		if obj.warnings
			warning("Disabling incomplete factorization, matrix is non-sparse.");
		end
	end
	
	% check solver method
	obj.method = arg.method;
	end
	%===========================================================
end

%% set/get methods
methods
	%===========================================================
	function set.method(obj,val)
	% update solver method
	arguments
		obj iterativeSolver
		val (1,1) string
	end

	% define valid method names
	validIterators = [ ...
		"lsqr","gmres","bicg","bicgstab","bicgstabl", ...
		"cgs","qmr","tfqmr","minres","symmlq","pcg" ...
		];
	
	if (val ~= obj.method) || ~matches(val,validIterators) || ~matches(obj.method,validIterators)
		% check the validity of the given solver method
		[obj.method] = iterativeSolver.check_method(obj.A,val,...
			"isSymm",obj.isSymm,"isHerm",obj.isHerm,"isPosDef",obj.isPosDef); %#ok

		% throw warning if the solver method was invalid
		if (val ~= obj.method) && (obj.warnings) %#ok
			warning("The suggested solver method for the iterativeSolver object " + ...
				"is being changed from """+val+ """ to """ + obj.method + """.");
		end
	end
	% DONE
	end
	%===========================================================
end
%% DONE
end