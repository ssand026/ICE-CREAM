function [varargout] = npermsk(n,k)
% For a set containing "n" distinct elements, return all possible permutations of length "k"
% The order of the elements of a permutation does matter: each row of the output
% contains "k" elements in the range [1,n], where repeated values are allowed.
%
% SEE ALSO: NCOMBSK
arguments
	n (1,1) {mustBeInteger,mustBePositive}
	k (1,1) {mustBeInteger,mustBePositive}
end

% try to initialize the permutation list
try 
	permList = zeros(n^k,k);
catch ME
	msg = "The requested permutation list is too large:"+newline+ME.message;
	throw(MException(ME.identifier,msg))
end

% populate the permutation list
for kk = 1:k
	permList(:,kk) = repmat(repelem((1:n).',n^(kk-1),1),n^(k-kk),1);
end

if nargout <= 1
	% output the list of permutations as an array
	varargout{1} = permList;
else
	% push each column of the permutation list to the output
	varargout = cell(1,nargout);
	for ii = 1:min(k,nargout)
		varargout{ii} = permList(:,ii);
	end
end
% DONE
end