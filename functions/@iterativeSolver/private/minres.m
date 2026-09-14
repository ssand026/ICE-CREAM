function [x,flag] = minres(A,b,tol,maxit,M1,M2,x0)
% MINRES Streamlined version of the Minimum Residual Method.
%
%   X = MINRES(A,B) attempts to find a minimum norm residual solution X to
%   the system of linear equations A*X=B. The N-by-N coefficient matrix A
%   must be symmetric but need not be positive definite. The right hand
%   side column vector B must have length N.
%
%   X = MINRES(A,B,TOL,MAXIT,M) and X = MINRES(A,B,TOL,MAXIT,M1,M2) use the
%   symmetric positive definite preconditioner M or M = M1*M2 and
%   effectively solve the symmetric system inv(H)*A*inv(H')*Y = inv(H)*B
%   for Y, where Y = H'*X and M = H*H'. The algorithm does not form H
%   explicitly. If M is [] then a preconditioner is not applied. M may be a
%   function handle MFUN returning M\X.
%
%   X = MINRES(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0
%   is [] then MINRES uses the default, an all zero vector.
%
%   [X,FLAG] = MINRES(A,B,...) also returns a convergence FLAG:
%    0 MINRES converged to the desired tolerance TOL within MAXIT
%    iterations.
%    1 MINRES iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 MINRES stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during MINRES became
%      too small or too large to continue computing.
%    5 preconditioner M was not symmetric positive definite.
%

%   Copyright 1984-2024 The MathWorks, Inc.
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
maxit = max(maxit,1);

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
stagSteps = 0;
maxStagSteps = 3;
moreSteps = 0;
maxMoreSteps = min([floor(numCol/50),5,numCol-maxit]);

% minimal residual solution
xMin = x;
rNormMin = rNorm;

% etc.
existM1 = ~isempty(M1);
existM2 = ~isempty(M2);
%-----------------------------------------------------------
v = r;
beta = 0;
vOld = v;
if existM1
	u = M1\vOld;
	if ~allfinite(u)
		flag = 2;
		return
	end
else % no preconditioner
	u = vOld;
end
if existM2
	v = M2\u;
	if ~allfinite(v)
		flag = 2;
		return
	end
else % no preconditioner
	v = u;
end
betaOld = beta;
beta = vOld' * v;
if (beta <= 0)
	flag = 5;
	return
end
beta = sqrt(beta);
%-----------------------------------------------------------
snprod = beta;
sn = 0;
cs = -1;
m = zeros(numCol,1);
mOld = m;
Am = zeros(numCol,1);
AmOld = Am;
deltaBar = beta;
epsilon = 0;

%% loop over maxit iterations (unless convergence or failure)
%===========================================================
for ii = 1:maxit
	%-----------------------------------------------------------
	vv = v / beta;
	v = A * vv;
	AmOlder = AmOld;
	AmOld = Am;
	Am = v;
	if ii ~= 1
		v = v - (beta / betaOld) * vOlder;
		alpha = vv' * v;
		v = v - (alpha / beta) * vOld;
	else
		alpha = vv' * v;
		v = v - (alpha/beta) * vOld;
		% Local reorthogonalization
		numer = vv' * v;
		denom = vv' * vv;
		v = v - (numer/denom) * vv;
	end
	%-----------------------------------------------------------
	vOlder = vOld;
	vOld = v;
	if existM1
		u = M1\vOld;
		if ~allfinite(u)
			flag = 2;
			return
		end
	else % no preconditioner
		u = vOld;
	end
	if existM2
		v = M2\u;
		if ~allfinite(v)
			flag = 2;
			return
		end
	else % no preconditioner
		v = u;
	end
	%-----------------------------------------------------------
	betaOld = beta;
	beta = vOld' * v;
	if (beta < 0)
		flag = 5;
		return
	end
	beta = sqrt(beta);
	%-----------------------------------------------------------
	delta = cs * deltaBar + sn * alpha;
	
	mOlder = mOld;
	mOld = m;
	m = vv - delta * mOld - epsilon * mOlder;
	Am = Am - delta * AmOld - epsilon * AmOlder;
	
	gammabar = sn * deltaBar - cs * alpha;
	epsilon = sn * beta;
	deltaBar = -cs * beta;
	gamma = sqrt(gammabar^2 + beta^2);
	
	m = m / gamma;
	Am = Am / gamma;
	cs = gammabar / gamma;
	sn = beta / gamma;
	%-----------------------------------------------------------
	% check for stagnation of the method
	if (snprod*cs == 0) || (abs(snprod*cs)*norm(m) < epsT*norm(x))
		% increment the number of consecutive iterates which are the same
		stagSteps = stagSteps + 1;
	else
		stagSteps = 0;
	end
	% This recurrence produces CG iterates.
	% Enable the following statement to see xcg.
	% xcg = x + snprod * (sn/cs) * m;
	x = x + snprod * cs * m;
	snprodOld = snprod;
	snprod = snprod * sn;
	if (existM1 || existM2)
		r = r - snprodOld * cs * Am;
		rNorm = norm(r);
	else
		rNorm = full(abs(snprod));
	end
	%-----------------------------------------------------------
	% check for convergence
	if (rNorm <= bTol) || (stagSteps >= maxStagSteps) || (moreSteps~=0)
		r = b - A*x;
		rNorm = norm(r);
		if rNorm <= bTol
			flag = 0;
			break
		else
			if (stagSteps >= maxStagSteps) && (moreSteps == 0)
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
	if (rNorm < rNormMin)
		rNormMin = rNorm;
		xMin = x;
	end
	% exit if stagnant
	if (stagSteps >= maxStagSteps)
		flag = 3;
		break
	end
	%-----------------------------------------------------------
end

%% return the minimal residual solution
%===========================================================
if (flag ~= 0)
	r_comp = b - A*xMin;
	if norm(r_comp) <= rNorm
		x = xMin;
	end
end
%% DONE
end