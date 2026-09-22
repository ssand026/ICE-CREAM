function [psi] = propagate_CN(psi, tau, H, M, opt)
% Forward-propagates a set of states using using the Crank-Nicolson propagator.
%
% SEE ALSO: ITERATIVESOLVER.LITSOLVE
arguments
	psi (:,:) double
	tau (1,1) double
	H double {mustBeSquare}
	M double {mustBeSquare} = [];
	% optional arguments
	opt.method (1,1) string = "builtin";
	opt.checkMethod (1,1) logical = false;
	% take iterativeSolver properties as additional arguments
	opt.tolerance (1,1) double {mustBeInRange(opt.tolerance,0,1)}
	opt.maxIters (1,1) uint16 {mustBePositive}
	opt.restarts (1,1) uint16
	opt.equilibrate (1,1) logical
	opt.factorize   (1,1) logical
	opt.useCholFact (1,1) logical
	opt.checkPosDef (1,1) logical
	opt.warnings (1,1) logical = false;
	opt.parallel (1,1) logical = false;
end

% check/parse inputs
%-----------------------------------------------------------
if isempty(M)
	% set overlap to the full identity matrix
	M = eye(size(H),"like",H);
end

% check that the inputs have appropriate sizes
if (size(H,2) ~= size(psi,1))
	error("ERROR: the sizes of 'H' and 'psi' must match")
elseif ~isequal(size(M),size(H))
	error("ERROR: the sizes of 'H' and 'M' must match")
end

% define valid methods for the iterative solver
useIterativeMethod = matches(opt.method,["iter",...
	"lsqr","gmres","bicg","bicgstab","bicgstabl", ...
	"cgs","qmr","tfqmr","minres","symmlq","pcg"]);

% apply the Crank-Nicolson operator
%-----------------------------------------------------------
if useIterativeMethod
	% use an iterative solver
	iterOpts = namedargs2cell(opt);
	tmp = (M - (1i*tau/2) * H) * psi;
	mat = (M + (1i*tau/2) * H);
	psi = iterativeSolver.litsolve(mat,tmp,psi,iterOpts{:});
else
	% use the default backslash operator
	psi = (M - (1i*tau/2) * H) * psi;
	psi = (M + (1i*tau/2) * H) \ psi;
end
% DONE
end