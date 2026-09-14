function [isSquare,isSymm,isHerm,isPosDef] = matrix_properties(A)
% determine key properties of the matrix 'A'
arguments
	A (:,:) double
end

% initialize property flags
isSquare = false;
isSymm = false;
isHerm = false;
isPosDef = false;

% check if matrix is square
if (nargout >= 1)
	isSquare = (size(A,1)==size(A,2));
end

% check if matrix is symmetric
if (nargout >= 2) && isSquare
	isSymm = issymmetric(A);
end
		
% check if matrix is hermitian
if (nargout >= 3) && isSquare
	if isSymm
		isHerm = isreal(A);
	else
		isHerm = ishermitian(A);
	end
end

% check if the matrix is positive-definite
if (nargout >= 4) && (isSymm || isHerm)
	try
		[~,cholFlag] = chol(A);
		isPosDef = (cholFlag == 0);
	catch ME
		isPosDef = false;
	end
end
% DONE
end