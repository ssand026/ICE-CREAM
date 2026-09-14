function [x,flag] = pcg(A,b,tol,maxit,M1,M2,x0)
% PCG Streamlined version of the Preconditioned Conjugate Gradients Method.
%
%   X = PCG(A,B) attempts to solve the system of linear equations A*X=B for
%   X. The N-by-N coefficient matrix A must be symmetric and positive
%   definite and the right hand side column vector B must have length N.
%
%   X = PCG(A,B,TOL,MAXIT,M) and X = PCG(A,B,TOL,MAXIT,M1,M2) use the 
%   symmetric positive definite preconditioner M or M = M1*M2 and
%   effectively solve the symmetric system inv(H)*A*inv(H')*Y = inv(H)*B
%   for Y, where Y = H'*X and M = H*H'. The algorithm does not form H
%   explicitly. If M is [] then a preconditioner is not applied. M may be a
%   function handle MFUN returning M\X.
%
%   X = PCG(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0 is
%   [] then PCG uses the default, an all zero vector.
%
%   [X,FLAG] = PCG(A,B,...) also returns a convergence FLAG:
%    0 PCG converged to the desired tolerance TOL within MAXIT iterations
%    1 PCG iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 PCG stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during PCG became too
%      small or too large to continue computing.
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

% other
existM1 = ~isempty(M1);
existM2 = ~isempty(M2);
rhoNew = 1;

%% loop over maxit iterations (unless convergence or failure)
%===========================================================
for ii = 1:maxit
	%-----------------------------------------------------------
	if existM1
		y = M1 \ r;
		if ~allfinite(y)
			flag = 2;
			break
		end
	else
		y = r;
	end
	if existM2
		z = M2 \ y;
		if ~allfinite(z)
			flag = 2;
			break
		end
	else
		z = y;
	end
	%-----------------------------------------------------------
	rhoOld = rhoNew;
	rhoNew = sum(r .* z,1);
	if (rhoNew == 0) || isinf(rhoNew)
		flag = 4;
		break
	end
	if (ii == 1)
		p = z;
		q = A * p;
	else
		beta = rhoNew / rhoOld;
		if (beta == 0) || isinf(beta)
			flag = 4;
			break
		end
		p = z + beta * p;
		q = A * p;
	end
	p_dot_q = p' * q;
	if (p_dot_q <= 0) || isinf(p_dot_q)
		flag = 4;
		break
	end
	alpha = rhoNew ./ p_dot_q;
	if isinf(alpha)
		flag = 4;
		break
	end
	%-----------------------------------------------------------
	% check for stagnation
	if (abs(alpha)*norm(p) < epsT*norm(x))
		stagSteps = stagSteps + 1;
	else
		stagSteps = 0;
	end
	% form new iterate
	x = x + alpha * p;
	r = r - alpha * q;
	rNorm = norm(r);
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
	if stagSteps >= maxStagSteps
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