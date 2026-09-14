function [x,flag] = solve(obj,b,x0)
% Solves A*x=b using the specified iterative solver. The argument 'x0'
% can be used to initialize the solution to the iterative solver if an
% approximate solution to the problem is known. If 'x0' is close to the
% ideal solution, this can reduce the time required to converge.
%
% The second output corresponds to the exit flags from the iterative
% solvers, which is typically zero if the solver converged.
%
% SEE ALSO:
%   EQUILIBRATE, ILU, ICHOL,
%   LSQR, GMRES, BICG, BICGSTAB, BICGSTABL, CGS, QMR, TFQMR, MINRES, SYMMLQ, PCG
arguments
	obj iterativeSolver
	b double {mustBeMatrix}
	x0 double {mustBeMatrix} = [];
end

% check input sizes
%-----------------------------------------------------------
[numRows,numCols] = size(obj.A);
numVec = size(b,2);
if size(b,1)~=numRows
	error("ERROR: the number of rows in the rhs vector '"+inputname(2)+ ...
		"' must equal the number of rows in the matrix '"+inputname(1)+"'.")
elseif ~isempty(x0) && (size(x0,1)~=numCols)
	error("ERROR: the number of rows in the initializtion vector '"+inputname(4) ...
		+"' must equal the number of columns in the matrix '"+inputname(1)+"'.")
elseif ~isempty(x0) && (size(x0,2)~=numVec)
	error("ERROR: the number of columns in the initializtion vector '"+inputname(4)+ ...
		"' must equal the number of columns in the rhs vector '"+inputname(2)+"'.")
end

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

% check for equillibrate preconditioners
%-----------------------------------------------------------
if obj.equilibrate
	if ~(isequal(size(obj.R),[numRows,1]) && isequal(size(obj.C),[numCols,1]))
		% valid equillibrate row/column scalings do not exist
		[rowPerm,obj.R,obj.C] = equilibrate(obj.A,"vector");
		[~,rowPerm] = sort(rowPerm);
		obj.R = obj.R(rowPerm,:);
	end
	
	if ~isequal(size(obj.B),[numRows,numCols])
		% a valid equilibrated matrix does not exist
		obj.B = obj.R .* obj.A .* obj.C.';
	end
	
	% apply preconditioners to the rhs/solution vectors
	b = obj.R .* b;
	if ~isempty(x0)
		x0 = x0 ./ obj.C;
	end

	% use equilibrated matrix going forward
	AA = obj.B;
else
	% use the non-equilibrated matrix going forward
	AA = obj.A;
end

% check if a valid factorization exists
%-----------------------------------------------------------
if ~obj.factorize
	% do not use an incomplete factorization
	M1 = [];
	M2 = [];
elseif obj.useCholFact && matches(obj.method,["pcg","minres","symmlq"])
	% use an incomplete Cholesky factorization
	if ~isequal(size(obj.Q),[numRows,numCols])
		% existing Cholesky factorization is invalid
		if obj.isPosDef
			% matrix is symmetric positive-definite, no diagcomp needed
			ichol_opts = struct("type",'nofill',"michol",'off');
			obj.Q = ichol(AA,ichol_opts);
		else
			% compute diagonal shift so that matrix is diagonally-dominant
			A_diag = abs(diag(AA));
			alpha = (sum(abs(AA),2) ./ A_diag) - 2;
			alpha = max(alpha(isfinite(alpha)));
			alpha = min(max(0,alpha),1/eps);
			ichol_opts = struct("type",'nofill',"michol",'off',"diagcomp",alpha);
			obj.Q = ichol(AA,ichol_opts);
		end
	end
	M1 = obj.Q;
	M2 = obj.Q';
	
else
	% use an incomplete LU factorization
	invalidLU = ...
		(size(obj.L,1)~=numRows) || ...
		(size(obj.U,2)~=numCols) || ...
		(size(obj.L,2)~=size(obj.U,1));
	if invalidLU
		% existing lu factorization is invalid
		ilu_opts = struct("type",'nofill',"milu",'off');
		[obj.L,obj.U] = ilu(AA,ilu_opts);
	end
	M1 = obj.L;
	M2 = obj.U;
end

% set solver method
%-----------------------------------------------------------
if (obj.method == "gmres")
	iterParams = {obj.restarts, obj.tolerance, obj.maxIters, M1, M2};
else
	iterParams = {obj.tolerance, obj.maxIters, M1, M2};
end

switch obj.method
	case "lsqr";      solvefun = @(b,x0)      lsqr(AA,b,iterParams{:},x0);
	case "gmres";     solvefun = @(b,x0)     gmres(AA,b,iterParams{:},x0);
	case "bicg";      solvefun = @(b,x0)      bicg(AA,b,iterParams{:},x0);
	case "bicgstab";  solvefun = @(b,x0)  bicgstab(AA,b,iterParams{:},x0);
	case "bicgstabl"; solvefun = @(b,x0) bicgstabl(AA,b,iterParams{:},x0);
	case "cgs";       solvefun = @(b,x0)       cgs(AA,b,iterParams{:},x0);
	case "qmr";       solvefun = @(b,x0)       qmr(AA,b,iterParams{:},x0);
	case "tfqmr";     solvefun = @(b,x0)     tfqmr(AA,b,iterParams{:},x0);
	case "minres";    solvefun = @(b,x0)    minres(AA,b,iterParams{:},x0);
	case "symmlq";    solvefun = @(b,x0)    symmlq(AA,b,iterParams{:},x0);
	case "pcg";       solvefun = @(b,x0)       pcg(AA,b,iterParams{:},x0);
end

% silence potential warnings
%-----------------------------------------------------------
if ~obj.warnings
	warning("off","MATLAB:"+obj.method+":tooSmallTolerance");
	warning("off","MATLAB:"+obj.method+":tooBigTolerance");
end

% solve the equation A*x = b for each column vector in b
%-----------------------------------------------------------
if isempty(x0); x0 = zeros(0,numVec); end
x = zeros(numCols,numVec);
flag = zeros(1,numVec);

if obj.parallel && ~isempty(gcp("nocreate"))
	parfor ii = 1:numVec
		[x(:,ii),flag(:,ii)] = solvefun(b(:,ii),x0(:,ii)); %#ok
	end
else
	% sequentially
	for ii = 1:numVec
		[x(:,ii),flag(:,ii)] = solvefun(b(:,ii),x0(:,ii));
	end
end

% undo equillibrate scaling
if obj.equilibrate
	x = obj.C .* x;
end

% restore potential warnings
%-----------------------------------------------------------
if ~obj.warnings
	warning("on","MATLAB:"+obj.method+":tooSmallTolerance");
	warning("on","MATLAB:"+obj.method+":tooBigTolerance");
end
% DONE
end