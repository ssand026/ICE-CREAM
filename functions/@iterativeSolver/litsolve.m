function [x,flag] = litsolve(A,b,x0,opt,pre)
% Solves A*x=b using the specified iterative solver. If the values of the
% solution vectors 'x' are approximately known, they can be input as the
% third argument 'x0' to speed up the convergence of the iterative solvers.
%
% The second output corresponds to the exit flags from the iterative
% solvers, which is typically zero if the solver converged.
%
% SEE ALSO:
%   ITERATIVESOLVER, EQUILIBRATE, ILU, ICHOL,
%   LSQR, GMRES, BICG, BICGSTAB, BICGSTABL, CGS, QMR, TFQMR, MINRES, SYMMLQ, PCG
arguments
	A  double {mustBeMatrix}
	b  double {mustBeMatrix}
	x0 double {mustBeMatrix} = [];
	% optional arguments
	opt.method (1,1) string
	opt.checkMethod (1,1) logical = false;
	% optional arguments in iterativeSolver
	opt.tolerance   (1,1) double {mustBeInRange(opt.tolerance,0,1)} = eps;
	opt.maxIters    (1,1) uint16 {mustBePositive} = 256;
	opt.restarts    (1,1) uint16 = 0;
	opt.equilibrate (1,1) logical = false;
	opt.factorize   (1,1) logical = true;
	opt.useCholFact (1,1) logical = false;
	opt.checkPosDef (1,1) logical = false;
	opt.warnings    (1,1) logical = true;
	opt.parallel    (1,1) logical = false;
	% preconditioner elements
	pre.R (:,1) double {mustBeVector} = []; % row-scaling vector from equilibrate
	pre.C (:,1) double {mustBeVector} = []; % column-scaling vector from equilibrate
	pre.B (:,:) double {mustBeMatrix} = []; % equilibrated matrix (R .* A .* C')
	pre.L (:,:) double {mustBeMatrix} = []; % lower ilu component
	pre.U (:,:) double {mustBeMatrix} = []; % upper ilu component
	pre.Q (:,:) double {mustBeMatrix} = []; % lower ichol component
end

% check matrix characteristics
isSparse = issparse(A);
[isSquare,isSymm,isHerm] = matrix_properties(A);

% check input sizes
%-----------------------------------------------------------
[numRows,numCols] = size(A);
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

% disable invalid preconditioners
%-----------------------------------------------------------
if (opt.equilibrate == true) && ~isSquare
	opt.equilibrate = false;
	if opt.warnings
		warning("Disabling equilibrate preconditioning, matrix is non-square.");
	end
end
if (opt.equilibrate == true) && (isSymm || isHerm)
	opt.equilibrate = false;
	if opt.warnings
		warning("Disabling equilibrate preconditioning, matrix is symmetric.");
	end
end
if (opt.factorize == true) && ~isSparse
	opt.factorize = false;
	if opt.warnings
		warning("Disabling incomplete factorization, matrix is non-sparse.");
	end
end


% check for equillibrate preconditioners
%-----------------------------------------------------------
if opt.equilibrate
	if ~(isequal(size(pre.R),[numRows,1]) && isequal(size(pre.C),[numCols,1]))
		% valid equillibrate row/column scalings do not exist
		[rowPerm,pre.R,pre.C] = equilibrate(A,"vector");
		[~,rowPerm] = sort(rowPerm);
		pre.R = pre.R(rowPerm,:);
	end
	
	if ~isequal(size(pre.B),[numRows,numCols])
		% a valid equilibrated matrix does not exist
		pre.B = pre.R .* A .* pre.C.';
	end
	
	% apply preconditioners to the rhs/solution vectors
	b = pre.R .* b;
	if ~isempty(x0)
		x0 = x0 ./ pre.C;
	end
	
	% use equilibrated matrix going forward
	AA = pre.B;
else
	% use the non-equilibrated matrix going forward
	AA = A;
end

% check solver method
%-----------------------------------------------------------
validIterators = [ ...
	"lsqr","gmres","bicg","bicgstab","bicgstabl", ...
	"cgs","qmr","tfqmr","minres","symmlq","pcg" ...
	];

% assume positive definiteness if matrix is symmetric/hermitian
isPosDef = (isSymm || isHerm);

if opt.checkMethod || ~matches(opt.method,validIterators)
	% check/update the given iterative method
	
	if opt.checkPosDef
		% actually check if matrix is positive-definite
		[~,~,~,isPosDef] = matrix_properties(AA);
	end
	
	% check the input method
	[method] = iterativeSolver.check_method(AA,opt.method,...
			"isSymm",isSymm,"isHerm",isHerm,"isPosDef",isPosDef);
	
	% throw warning if the solver method was invalid
	if opt.warnings && (method ~= opt.method)
		warning("The suggested solver method for the itsolver object " + ...
			"is being changed from """+opt.method+ """ to """ + method + """.");
	end
else
	% use the given solver method
	method = opt.method;
end


% check if a valid factorization exists
%-----------------------------------------------------------
if ~opt.factorize
	% do not use an incomplete factorization
	M1 = [];
	M2 = [];
elseif opt.useCholFact && matches(method,["pcg","minres","symmlq"])
	% use an incomplete Cholesky factorization
	if ~isequal(size(pre.Q),[numRows,numCols])
		% existing Cholesky factorization is invalid
		if isPosDef
			% matrix is symmetric positive-definite, no diagcomp needed
			ichol_opts = struct("type",'nofill',"michol",'off');
			pre.Q = ichol(AA,ichol_opts);
		else
			% compute diagonal shift so that matrix is diagonally-dominant
			A_diag = abs(diag(AA));
			alpha = (sum(abs(AA),2) ./ A_diag) - 2;
			alpha = max(alpha(isfinite(alpha)));
			alpha = min(max(0,alpha),1/eps);
			ichol_opts = struct("type",'nofill',"michol",'off',"diagcomp",alpha);
			pre.Q = ichol(AA,ichol_opts);
		end
	end
	M1 = pre.Q;
	M2 = pre.Q';
	
else
	% use an incomplete LU factorization
	invalidLU = ...
		(size(pre.L,1)~=numRows) || ...
		(size(pre.U,2)~=numCols) || ...
		(size(pre.L,2)~=size(pre.U,1));
	if invalidLU
		% existing lu factorization is invalid
		ilu_opts = struct("type",'nofill',"milu",'off');
		[pre.L,pre.U] = ilu(AA,ilu_opts);
	end
	M1 = pre.L;
	M2 = pre.U;
end

% set solver method
%-----------------------------------------------------------
if (method == "gmres")
	iterParams = {opt.restarts, opt.tolerance, opt.maxIters, M1, M2};
else
	iterParams = {opt.tolerance, opt.maxIters, M1, M2};
end

switch method
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
	otherwise
		error("ERROR: invalid iterative solver name");
end

% silence potential warnings
%-----------------------------------------------------------
if ~opt.warnings
	warning("off","MATLAB:"+method+":tooSmallTolerance");
	warning("off","MATLAB:"+method+":tooBigTolerance");
end

% solve the equation A*x = b for each column vector in b
%-----------------------------------------------------------
if isempty(x0); x0 = zeros(0,numVec); end
x = zeros(numCols,numVec);
flag = zeros(1,numVec);

if opt.parallel && ~isempty(gcp("nocreate"))
	% in parallel
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
if opt.equilibrate
	x = pre.C .* x;
end

% restore potential warnings
if ~opt.warnings
	warning("on","MATLAB:"+method+":tooSmallTolerance");
	warning("on","MATLAB:"+method+":tooBigTolerance");
end
% DONE
end