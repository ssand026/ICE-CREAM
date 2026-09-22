function [H] = drop_perturb(H, drop_tol)
% DROP_PERTURB Drops minor contributions from a matrix using an energy-based criterion
% 
% Returns a dropped form of the matrix H, where the element H_ij -> 0 if the
% relative contribution to the total energy is less than "drop_tol"
%
% Useful for working in an orthonormal basis as this approximation can restore the sparsity of the
% Hamiltonian/operators of a system, reducing compute costs.
%
% SEE ALSO: SPARSIFY
arguments
	H double {mustBeSquare}
	drop_tol (1,1) double {mustBeInRange(drop_tol,0,1)}
end
if issparse(H); H = full(H); end
if drop_tol==0; return; end

% get the energy cutoff from the spectral range of the eigenvalues
Hjj = diag(H);
cutoff = drop_tol * abs(max(Hjj) - min(Hjj));
energy_contrib = abs(H .* H') ./ abs(Hjj.' - Hjj);

% drop values below the cutoff
to_drop = (energy_contrib < cutoff);
H(to_drop) = 0;
end