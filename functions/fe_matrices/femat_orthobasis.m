function [intoB,exitB] = femat_orthobasis(mesh,BC)
% Returns the orthonormal basis generators for a finite-element system with the
% given boundary conditions.
%
% The output 'intoB' is the basis reduction matrix
% The output 'exitB' is the basis expansion matrix
%
% Calculations done in the reduced basis will automatically satisfy the
% given boundary conditions while having a identity overlap-matrix.
arguments
	mesh FE_mesh
	BC double {mustBeMatrix} = femat_dirichlet(mesh);
end

% apply the boundary conditions to the overlap matrix
M = FEmat(mesh,"overlap");
if ~isempty(BC)
	M = to_basis(M,"reduced",BC);
end

% get the orthonormal-basis reversal matrix
exitB = sqrtm(full(M));

% get the orthonormal-basis application matrix
intoB = inv(exitB);
intoB = (intoB + intoB.')/2;

% apply the basis expansion/reduction to the matrix
if ~isempty(BC)
	if (size(BC,1) > size(BC,2)); BC = BC.'; end
	intoB = intoB * BC;
	exitB = BC.' * exitB;
end
% DONE
end