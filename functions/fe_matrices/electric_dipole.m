function [muE] = electric_dipole(mesh,perm,q)
% Constructs the dipole-moment matrices for spatially-invariant electric fields.
% The output will be a cell array, where each entry corresponds to the dipole
% operator for the field along each axis.
arguments
	mesh FE_mesh
	perm {mustBeA(perm,{'FE_coeff','function_handle','char','string','numeric'})} = 1;
	q (1,1) double {mustBeNumeric,mustBeFinite} = -1;
end

% convert permitivity to a coefficient
if ~isa(perm,'FE_coeff')
	perm = FE_coeff(mesh,perm);
end

% define the position vector
r = FE_coeff(mesh,"@(pts) pts");
r = r.';

% compute electric dipole moments
[numDim] = mesh.num("dim");
muE = cell(1,numDim);
E_field = eye(numDim,numDim);

for ii = 1:numDim
	D_ii = perm * E_field(:,ii);
	V_ii = r * D_ii;
	muE{ii} = q * FEmat(mesh,"scalar",V_ii);
end
% DONE
end