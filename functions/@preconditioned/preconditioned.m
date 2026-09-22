classdef preconditioned < handle & matlab.mixin.Copyable & matlab.mixin.SetGet & ...
		matlab.mixin.CustomElementSerialization
	% Preconditions a matrix to enable efficient repeated solves of a system of linear
	% equations, as it eliminates repeated calculation of any preconditioning steps.
	% While this can often amount to using the builtin "decomposition" function, for
	% small dense matrices, direct matrix inversion may be more efficient, while for
	% large sparse matrices, iterative methods will be better
	%
	% SOLVER METHODS
	%  - "invert": explicit matrix inversion, only useful for small/dense matrices
	%  - "decompose": uses the builtin decomposition function, generally the best (non-iterative) approach
	%  - "iter": will automatically choose the iterative solver method
	%  - "lsqr": only iterative method that works for non-square matrices
	%  - "gmres": is typically the most reliable/flexible iterative solver method
	%  - "minres": good for symmetric matrices
	%  - "cgs":  can be faster than minres/gmres in certain cases
	%  - "pcg": good for symmetric positive-definite matrices
	%
	%  Other iterative solvers have fairly specific use cases
	%
	% SEE ALSO: DECOMPOSITION, ITERATIVESOLVER, EQUILIBRATE, ILU, ICHOL,
	%   LSQR, GMRES, BICG, BICGSTAB, BICGSTABL, CGS, QMR, TFQMR, MINRES, SYMMLQ, PCG
	
	
	%% properties
	properties (SetAccess=protected)
		A double {mustBeMatrix}
		isSparse (1,1) logical
		isSquare (1,1) logical
		isSymm (1,1) logical
		isHerm (1,1) logical
		isPosDef (1,1) logical
	end
	
	properties (SetObservable)
		method (1,1) string = "decompose";
	end
	
	% iterativeSolver options
	properties
		tolerance (1,1) double {mustBeInRange(tolerance,0,1)} = 1e-12;
		maxIters (1,1) uint16 {mustBePositive,mustBeInteger} = 256;
		restarts (1,1) uint16 {mustBeInteger} = 0;
		equilibrate (1,1) logical = false;
		factorize   (1,1) logical = true;
		useCholFact (1,1) logical = false;
		checkPosDef (1,1) logical = false;
		warnings (1,1) logical = false;
		parallel (1,1) logical = false;
	end
	
	% preconditioner objects
	properties (Hidden)
		dA = []
		iA = []
		lA = []
		rA = []
	end

	% prevent saving/loading issues with the decomposition class
	methods (Static,Hidden)
		%===========================================================
		function modifyOutgoingSerializationContent(saveObj,~,~)
		% clear decomposition objects before saving
		saveObj.dA = [];
		end
		%===========================================================
	end
	
	%% constructor method
	methods
		%===========================================================
		function [obj] = preconditioned(A,arg,opt)
		arguments
			A double {mustBeMatrix}
			arg.method (1,1) {mustBeMember(arg.method,["decompose","invert","iter", ...
				"lsqr","gmres","bicg","bicgstab","bicgstabl","cgs","qmr","tfqmr","minres","symmlq","pcg"])}
			% optional properties for iterativeSolver
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
		
		% update matrix property values
		obj.A = A;
		obj.isSparse = issparse(A);
		obj.isSquare = isSquare;
		obj.isSymm = isSymm;
		obj.isHerm = isHerm;
		obj.isPosDef = isPosDef;
		
		% check solver method
		obj.method = arg.method;
		
		% evaluate pre-conditioners for the current solver method
		if (obj.method=="decompose")
			obj.dA = decomposition(A);
		elseif (obj.method == "invert")
			obj.iA = inv(A);
		else
			iterOpts = { ...
				"tolerance",obj.tolerance,"maxIters",obj.maxIters, ...
				"restarts",obj.restarts,"equilibrate",obj.equilibrate, ...
				"factorize",obj.factorize,"useCholFact",obj.useCholFact, ...
				"checkPosDef",obj.checkPosDef,"warnings",obj.warnings, ...
				"parallel",obj.parallel,"method",obj.method ...
				};
			obj.lA = iterativeSolver(A, iterOpts{:});
			if obj.isSymm
				obj.rA = obj.lA;
			else
				obj.rA = iterativeSolver(A.',iterOpts{:});
			end
		end
		% DONE
		end
		%===========================================================
	end
	
	%% set/get methods
	methods
		%===========================================================
		function set.method(obj,val)
		% update solver method
		arguments
			obj
			val (1,1) {mustBeMember(val,["decompose","invert","iter", ...
				"lsqr","gmres","bicg","bicgstab","bicgstabl","cgs","qmr","tfqmr","minres","symmlq","pcg"])}
		end
		
		if (val == obj.method)
			% input is already a valid method, do nothing
		elseif matches(val,["decompose","invert"])
			if ~obj.isSquare %#ok
				% switch to decomposition for non-square matrices
				obj.method = "decompose";
			else
				% use inputted method
				obj.method = val;
			end
		else
			% check the iterative solver method
			[obj.method] = iterativeSolver.check_method(obj.A,val,...
				"isSymm",obj.isSymm,"isHerm",obj.isHerm,"isPosDef",obj.isPosDef); %#ok
		end
		
		% throw warning if the solver method was invalid
		if (val ~= obj.method) && (obj.warnings) %#ok
			warning("The suggested solver method for the itsolver object " + ...
				"is being changed from """+val+ """ to """ + obj.method + """.");
		end
		% DONE
		end
		%===========================================================
	end
	
	%% operator overload
	methods
		%===========================================================
		function [out] = size(obj)
		out = size(obj.A);
		end
		%===========================================================
		function [out] = double(obj)
		out = double(obj.A);
		end
		%===========================================================
		function [out] = mtimes(A,B)
		if isa(A,"preconditioned"); A = double(A); end
		if isa(B,"preconditioned"); B = double(B); end
		out = A * B;
		end
		%===========================================================
		function [x] = mldivide(obj,b)
		arguments
			obj preconditioned
			b double
		end
		% perform linear solve using ...
		switch obj.method
			case "decompose"
				% the matrix decomposition
				if isempty(obj.dA)
					obj.dA = decomposition(obj.A);
				end
				x = obj.dA \ b;
			case "invert"
				% the matrix inverse
				if isempty(obj.iA)
					obj.iA = inv(obj.A);
				end
				x = obj.iA * b;
			otherwise
				% an iterative solver
				iterOpts = { ...
					"method",obj.method, ...
					"tolerance", obj.tolerance, ...
					"maxIters", obj.maxIters, ...
					"restarts", obj.restarts, ...
					"equilibrate", obj.equilibrate, ...
					"factorize", obj.factorize, ...
					"useCholFact", obj.useCholFact, ...
					"checkPosDef", obj.checkPosDef, ...
					"warnings", obj.warnings, ...
					"parallel", obj.parallel ...
					};
				if isempty(obj.lA)
					obj.lA = iterativeSolver(obj.A,iterOpts{:});
				else
					set(obj.lA,iterOpts(1:2:end),iterOpts(2:2:end));
				end
				x = solve(obj.lA,b);
		end
		% DONE
		end
		%===========================================================
		function [x] = mrdivide(b,obj)
		arguments
			b double {mustBeMatrix}
			obj preconditioned
		end
		% perform linear solve using ...
		switch obj.method
			case "decompose"
				% the matrix decomposition
				if isempty(obj.dA)
					obj.dA = decomposition(obj.A);
				end
				x = b / obj.dA;
			case "invert"
				% the matrix inverse
				if isempty(obj.iA)
					obj.iA = inv(obj.A);
				end
				x = b * obj.iA;
			otherwise
				% an iterative solver
				iterOpts = { ...
					"method",obj.method, ...
					"tolerance", obj.tolerance, ...
					"maxIters", obj.maxIters, ...
					"restarts", obj.restarts, ...
					"equilibrate", obj.equilibrate, ...
					"factorize", obj.factorize, ...
					"useCholFact", obj.useCholFact, ...
					"checkPosDef", obj.checkPosDef, ...
					"warnings", obj.warnings, ...
					"parallel", obj.parallel ...
					};
				if isempty(obj.rA)
					obj.rA = iterativeSolver(obj.A.',iterOpts{:});
				else
					set(obj.rA,iterOpts(1:2:end),iterOpts(2:2:end));
				end
				x = solve(obj.rA,b).';
		end
		% DONE
		end
		%===========================================================
	end
	%% DONE
end
