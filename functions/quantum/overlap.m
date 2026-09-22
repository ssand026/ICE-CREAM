function [P] = overlap(bra,M,ket)
%OVERLAP: Computes the overlap between a set of bra and ket vectors
arguments
	bra (:,:,:) double {mustBeNonsparse}
	M double {mustBeSquare}
	ket (:,:,:) double {mustBeNonsparse}
end

% check bra/ket sizes
[numBra,numRow] = size(bra);
[numCol,numKet] = size(ket);

minPages = min(size(bra,3),size(ket,3));
if (size(bra,3)~=size(ket,3)) && (minPages~=1)
	error("ERROR: the third dimension of bra and ket have incompatible lengths.")
end

% check if the overlap matrix is the identity matrix
[~,scale] = iseye(M);
isEye = ~isnan(scale);

if isEye
	% orthonormal basis
	%-----------------------------
	if (numRow ~= numCol)
		error("ERROR: the sizes of bra and ket do not agree.")
	end

	% compute the expectation values <bra|M|ket>
	sqrtP = scale * pagemtimes(bra,ket);
else
	% non-orthonormal basis
	%-----------------------------
	if (size(M,1) ~= numRow) || (size(M,2) ~= numCol)
		error("ERROR: the sizes of the bra and ket do not match the overlap matrix.")
	end
	
	% compute the expectation values <bra|M|ket>
	if ~issparse(M)
		sqrtP = pagemtimes(bra,pagemtimes(M,ket));
	elseif (minPages > 1)
		sqrtP = zeros(numBra,numKet,minPages);
		for ii = 1:minPages
			sqrtP(:,:,ii) = bra(:,:,ii) * M * ket(:,:,ii);
		end	
	elseif ismatrix(bra)
		sqrtP = pagemtimes(bra*M,ket);
	elseif ismatrix(ket)
		sqrtP = pagemtimes(bra,M*ket);
	end
end

% return |<bra|ket>|^2 and ensure P is bound between [0,1]
%----------------------------------------------------------
P = abs(sqrtP).^2;
P = max(0,min(P,1));

% DONE
end