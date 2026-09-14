function [isEye,scale] = iseye(X)
% ISEYE(X) check if the input is an identity matrix
%   
% OUTPUTS
%   isEye: equals true if X is an identity matrix, false otherwise
%   scale: returns the identity matrix scaling, is set to NaN if X 
%          is not a contant-valued square diagonal matrix 
%
isEye = false;
scale = NaN;
if (isnumeric(X) || islogical(X)) && ismatrix(X) && (size(X,1)==size(X,2)) && isdiag(X)
	% input is a square, diagonal, matrix
	if isempty(X)
		scale = 1;
		isEye = true;
	elseif all(diag(X)==X(1,1))
		% diagonals are same-valued, check diagonal values
		scale = double(X(1,1));
		isEye = isequal(scale,1);
	end
end
% DONE
end