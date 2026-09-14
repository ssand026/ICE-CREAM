function [newMethod] = check_method(A,method,opt)
% Returns/validates a iterative solver method for the matrix 'A'. 
% 
% The third input 'isPosDef' denotes if the matrix 'A' is positive definite. If
% 'isPosDef' is not specified or left empty, the function will try to compute 
% the Cholesky factorization of 'A' to determine if it is positive definite.
% Setting this value explictly to true/false will skip the expensive Cholesky
% factorization step for square and symmetric matrices.
%
% SEE ALSO: 
%   LSQR, GMRES, BICG, BICGSTAB, BICGSTABL, CGS, QMR, TFQMR, MINRES, SYMMLQ, PCG, CHOL
arguments
	A (:,:) double
	method (1,1) string = "unset";
	opt.isSymm (1,:) logical {mustBeScalarOrEmpty} = [];
	opt.isHerm (1,:) logical {mustBeScalarOrEmpty} = [];
	opt.isPosDef (1,:) logical {mustBeScalarOrEmpty} = [];
end

% check the method category
%-----------------------------------------------------------
% the given method name is valid when the matrix is ...
if matches(method,["lsqr"])
	% any matrix
	mtype = 1;
elseif matches(method,["gmres","bicg","bicgstab","bicgstabl","cgs","qmr","tfqmr"])
	% square
	mtype = 2;
elseif matches(method,["minres","symmlq"])
	% square, symmetric/hermitian
	mtype = 3;
elseif matches(method,["pcg"])
	% square, symmetric/hermitian, positive-definite
	mtype = 4;
else
	% the method name is invalid
	mtype = 5;
end

% check the properties of the matrix
%-----------------------------------------------------------
if (mtype > 1)
	% check if matrix is square
	isSquare = (size(A,1)==size(A,2));
end

if (mtype > 2)
 	% check if matrix is symmetric
	if ~isSquare 
		isSymm = false;
	elseif ~isempty(opt.isSymm)
		isSymm = opt.isSymm;
	else
		isSymm = issymmetric(A);
	end

	% check if matrix is hermitian
	if ~isSquare
		isHerm = false;
	elseif ~isempty(opt.isHerm)
		isHerm = opt.isHerm;
	elseif isSymm
		isHerm = isreal(A);
	else
		isHerm = ishermitian(A);
	end
end

if (mtype > 3)
	% check if matrix is positive definite
	if ~isSquare || (~isSymm && ~isHerm)
		isPosDef = false;
	elseif ~isempty(opt.isPosDef)
		isPosDef = opt.isPosDef;
	else
		try
			[~,cholFlag] = chol(A);
			isPosDef = (cholFlag == 0);
		catch ME
			isPosDef = false;
		end
	end
end

% switch to a valid solver method if the given method is invalid
%-----------------------------------------------------------
if (mtype > 1) && ~isSquare
	% switch to general solver
	newMethod = "lsqr";
elseif (mtype > 2) && ~(isSymm || isHerm)
	% switch to square solver
	newMethod = "gmres";
elseif (mtype > 3) && ~isPosDef
	% switch to symmetric/hermitian solver
	newMethod = "minres";
elseif (mtype > 4) && isPosDef
	% switch to a positive-definite solver
	newMethod = "pcg";
else
	% the original method was valid
	newMethod = method;
end
% DONE
end


