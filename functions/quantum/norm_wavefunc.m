function [psi,scl] = norm_wavefunc(psi,M,dim)
% normalizes the given set of wavefunctions 'psi', where the state-vectors lie
% along dimension 'dim'.
arguments
	psi double
	M double {mustBeSquare} = [];
	dim (1,1) {mustBeInteger,mustBePositive} = 1;
end

if iseye(M)
	% orthonormal basis
	%-----------------------------
	scl = vecnorm(abs(psi),2,dim);
	psi = psi ./ scl;
else
	% non-orthonormal basis
	%-----------------------------
	if (length(M) ~= size(psi,dim))
		error("ERROR: the sizes of the wavefunction and the overlap matrix must agree.")
	end
	
	% bring contracted dimension to the front
	dimOrder = [dim,1:dim-1,dim+1:ndims(psi)];
	sz = size(psi,dimOrder);
	psi = reshape(permute(psi,dimOrder),sz(1),prod([1,sz(2:end)]));
	
	% apply the normalization scaling factor
	scl = sqrt(sum( abs(conj(psi) .* (M * psi)), 1));
	psi = psi ./ scl;
	
	% restore original size/shape
	scl = ipermute(reshape(scl,[1,sz(2:end)]),dimOrder);
	psi = ipermute(reshape(psi,sz),dimOrder);
end
% DONE
end