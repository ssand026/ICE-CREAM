function [out] = near_hermitian(A,tol)
% NEAR_HERMITIAN Determine if a matrix is hermitian within the given tolerance.
%
% Input Arguments
% ================
%	A: input matrix
%	tol (optional): value between 0 and 1 
% 
% Output Arguments
% =================
%	+1 if A is hermitian
%	-1 if A is skew-hermitian
%	 0 if A is non-hermitian
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
isHerm = ishermitian(A,"nonskew");
isSkew = ishermitian(A,"skew");

if (tol==0) && (~isHerm && ~isSkew)
	% input does not have exact symmetry
	out = 0;
	return
elseif isHerm
	% input is exactly hermitian
	out = +1;
	return
elseif isSkew
	% input is exactly skew-hermitian
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

% check if input is hermitian
if ~isHerm && ~isSkew
	A_herm = abs(A - A');
	if all(A_herm < min_err,"all")
		isHerm = true;
	elseif any(A_herm > max_err,"all")
		isHerm = false;
	else
		isHerm = all(A_herm < error,"all");
	end
end

% check if input is skew-hermitian
if ~isHerm && ~isSkew
	A_skew = abs(A + A');
	if all(A_herm < min_err,"all")
		isSkew = true;
	elseif any(A_skew > max_err,"all")
		isSkew = false;
	else
		isSkew = all(A_skew < error,"all");
	end
end

if ~isHerm && ~isSkew
	% input is not hermitian
	out = 0;
elseif isHerm
	% input is hermitian
	out = +1;
elseif isSkew
	% input is skew-hermitian
	out = -1;
end
end