function [varargout] = ncombsk(n,k,opt)
% For a set containing "n" distinct elements, return all possible combinations of length "k"
% The order of the elements in a combination does NOT matter: each row of the output
% contains "k" elements in the range [1,n] without repeats.
%
% SEE ALSO: NPERMSK
arguments
	n (1,1) uint8
	k (1,1) uint8
	opt.asIndex (1,1) logical = false;
end

% get the logical index for elements in the combination
if     k==0;  combIndex = false(1,n);
elseif k==1;  combIndex = logical(eye(n));
elseif k==n;  combIndex = true(1,n);
else
	% evaluate combIndex recursively
	vals = cell(n,k);
	vals = cellfun(@logical,vals,'UniformOutput',false);
	for nn = 1:n
		for kk = 1:min(nn,k)
			if     kk==0;  vals{nn,kk} = false(1,nn);
			elseif kk==1;  vals{nn,kk} = logical(eye(nn));
			elseif kk==n;  vals{nn,kk} = true(1,nn);
			else
				% combine previous values of combIndex
				A = vals{nn-1,kk-1}; numA = size(A,1);
				B = vals{nn-1,kk+0}; numB = size(B,1);
				vals{nn,kk} = [[true(numA,1),A]; [false(numB,1),B]];
			end
		end
	end
	combIndex = vals{n,k};
end

% choose output form
if opt.asIndex
	% return the indices in logical form
	combList = logical(combIndex);
else
	% convert the logical index to a list of elements
	[combList,~] = find(permute(combIndex,[2,1]));
	combList = permute(reshape(combList,k,[]),[2,1]);
end

% prepare output
if nargout<=1
	% return the full list of combinations
	varargout{1} = combList;
else
	% push each column of the combination list to the output
	varargout = cell(1,nargout);
	for ii = 1:min(k,nargout)
		varargout{ii} = combList(:,ii);
	end
end
% DONE
end