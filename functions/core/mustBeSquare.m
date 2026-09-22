function mustBeSquare(A)
% MUSTBESQUARE Validate that the value is a square matrix
%	MUSTBESQUARE(A) throws an error if A is not a matrix or is not square.
%	empty. Calls ismatrix to determine if A is a matrix, and calls size
%	to determine if A is square.

if ~isnumeric(A) || ~ismatrix(A) || (size(A,1)~=size(A,2))
	ME = MException("mustBeSquare:validation","Input must be a square matrix");
	throwAsCaller(ME);
end
% DONE
end
