function [out] = near_symmetric(A,tol)
% NEAR_SYMMETRIC Determine if a matrix is symmetric within the given tolerance.
%
% Input Arguments
% ================
%	A: input matrix
%	tol (optional): value between 0 and 1 
% 
% Output Arguments
% =================
%	+1 if A is symmetric
%	-1 if A is skew-symmetric
%	 0 if A is non-symmetric
%
arguments
	A   double {mustBeMatrix}
	tol (1,1) {mustBeInRange(tol,0,1)} = 1e-12;
end

if size(A,1)~=size(A,2)
	% input matrix is not square
	out = 0;
	return;
end

% check if input has exact symmetry
isSymm = issymmetric(A,"nonskew");
isSkew = issymmetric(A,"skew");

if (~isSymm && ~isSkew)
	% input is not exactly symmetric
	if tol==0
		out = 0;
		return
	end
elseif isSymm
	% input is exactly symmetric
	out = +1;
	return
elseif isSkew
	% input is exactly skew-symmetric
	out = -1;
	return
end

% find the error threshold for the input
if issparse(A)
	error = (tol/eps)*spfun(@eps,A);
	error = (error + error.');
	max_err = max(nonzeros(error));
	min_err = min(nonzeros(error));
else
	error = (tol/eps)*eps(A);
	error = (error + error.');
	max_err = max(error(:));
	min_err = min(error(:));
end

% check if input is symmetric
if ~isSymm && ~isSkew
	A_symm = abs(A - A.');

	if all(A_symm < min_err,"all")
		isSymm = true;
	elseif any(A_symm > max_err,"all")
		isSymm = false;
	else
		isSymm = all(A_symm < error,"all");
	end
end

% check if input is skew-symmetric
if ~isSymm && ~isSkew
	A_skew = abs(A + A.');

	if all(A_symm < min_err,"all")
		isSkew = true;
	elseif any(A_skew > max_err,"all")
		isSkew = false;
	else
		isSkew = all(A_skew < error,"all");
	end
end

if ~isSymm && ~isSkew
	% input is not symmetric
	out = 0;
elseif isSymm
	% input is symmetric
	out = +1;
elseif isSkew
	% input is skew-symmetric
	out = -1;
end
end