function [x,flag] = cgs(A,b,tol,maxit,M1,M2,x0)
% CGS Streamlined version of the Conjugate Gradients Squared Method.
%
%   X = CGS(A,B) attempts to solve the system of linear equations A*X=B for
%   X. The N-by-N coefficient matrix A must be square and the right hand
%   side column vector B must have length N.
%
%   X = CGS(A,B,TOL,MAXIT,M) and X = CGS(A,B,TOL,MAXIT,M1,M2) use the
%   preconditioner M or M=M1*M2 and effectively solve the system
%   A*inv(M)*Y = B for Y, where Y = M*X. If M is [] then a preconditioner
%   is not applied.  M may be a function handle returning M\X.
%
%   X = CGS(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0 is
%   [] then CGS uses the default, an all zero vector.
%
%   [X,FLAG] = CGS(A,B,...) also returns a convergence FLAG:
%    0 CGS converged to the desired tolerance TOL within MAXIT iterations.
%    1 CGS iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 CGS stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during CGS became too
%      small or too large to continue computing.
%
%   See also BICG, BICGSTAB, BICGSTABL, GMRES, LSQR, MINRES, PCG, QMR,
%   SYMMLQ, TFQMR, ILU, FUNCTION_HANDLE.

%   Copyright 1984-2025 The MathWorks, Inc.
arguments
	A (:,:) double {mustBeSquare}
	b (:,:) double
	tol   (1,1) double = 1e-12;
	maxit (1,1) uint16 = 32;
	M1 (:,:) double = [];
	M2 (:,:) double = [];
	x0 (:,:) double = [];
end

%% validate inputs/intialize outputs
%===========================================================
% check input sizes
[numRow,numCol] = size(A);
if (size(b,1)~=numRow)
	error("ERROR: the number of rows in the right-hand-side vector 'b' must " + ...
		"match the number of columns in the matrix 'A'.");
end

% ensure input values are valid
epsT = eps("like",b);
tol = max(tol,epsT);
tol = min(tol,1-epsT);
maxit = min(maxit,numCol);

% initialize the solution vector
if isempty(x0)
	x = zeros([numCol,1],"like",b);
elseif isequal(size(x0),[numCol,1])
	x = x0;
else
	% x0 has invalid size
	error("ERROR: the size of the intial solution vector 'x0' is incompatible " + ...
		"with the inputs 'A' and 'b'.");
end

% initialize the convergence flag
flag = 1;

%% check for possible exit conditions
%===========================================================
% check if right-hand side is all zeros
bNorm = norm(b);
bTol = tol * bNorm;
if (bNorm==0)
	flag = 0;
	x(:) = 0;
	return
end

% check if initial guess is within required tolerance
r = b - A*x; % residual vector
rNorm = norm(r);
if (rNorm <= bTol)
	flag = 0;
	return
end

%% initialize loop variables
%===========================================================
% exit parameters
stagSteps = 0; % number of stagnant iterations
maxStagSteps = 3;
moreSteps = 0;
maxMoreSteps = min([floor(numCol/50),5,numCol-maxit]);

% minimal residual solution
xmin = x;
rNormMin = rNorm;

% other
existM1 = ~isempty(M1);
existM2 = ~isempty(M2);
rhoNew = 1;
rOrig = r; % original residual
u = r;
p = r;

%% loop over maxit iterations (unless convergence or failure)
%===========================================================
for ii = 1:maxit
	%-----------------------------------------------------------
	rhoOld = rhoNew;
	rhoNew = rOrig' * r;
	if (rhoNew == 0) || isinf(rhoNew)
		flag = 4;
		break
	end
	if (ii ~= 1)
		beta = rhoNew / rhoOld;
		if (beta == 0) || isinf(beta)
			flag = 4;
			break
		end
		u = r + beta * q;
		p = u + beta * (q + beta * p);
	end
	%-----------------------------------------------------------
	ph = p;
	if existM1
		ph = M1 \ ph;
	end
	if existM2
		ph = M2 \ ph;
	end
	if (existM1 || existM2) && ~allfinite(ph)
		flag = 2;
		break
	end
	ph = A * ph;
	%-----------------------------------------------------------
	alpha = rhoNew / (rOrig' * ph);
	if ~isfinite(alpha)
		flag = 4;
		break
	end
	q = u - alpha * ph;
	uh = u + q;
	if existM1
		uh = M1 \ uh;
	end
	if existM2
		uh = M2 \ uh;
	end
	if (existM1 || existM2) && ~allfinite(uh)
		flag = 2;
		break
	end
	%-----------------------------------------------------------
	% check for stagnation
	if abs(alpha)*norm(uh) < epsT*norm(x)
		stagSteps = stagSteps + 1;
	else
		stagSteps = 0;
	end
	% form the new iterate
	rh = (A * uh);
	x = x + alpha * uh;
	r = r - alpha * rh;
	rNorm = norm(r);
	%-----------------------------------------------------------
	% check for convergence
	if (rNorm <= bTol || stagSteps >= maxStagSteps || moreSteps)
		r = b - A * x;
		rNorm = norm(r);
		if (rNorm <= bTol)
			flag = 0;
			break
		else
			if stagSteps >= maxStagSteps && moreSteps == 0
				stagSteps = 0;
			end
			moreSteps = moreSteps + 1;
			if moreSteps >= maxMoreSteps
				flag = 3;
				break
			end
		end
	end
	% update minimal norm quantities
	if rNorm < rNormMin
		rNormMin = rNorm;
		xmin = x;
	end
	% exit if stagnant
	if stagSteps >= maxStagSteps
		flag = 3;
		break
	end
end

%% return the minimal residual solution
%===========================================================
if flag ~= 0
	r = b - A * xmin;
	if norm(r) <= rNorm
		x = xmin;
	end
end
%% DONE
end