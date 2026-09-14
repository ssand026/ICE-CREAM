function [muE,muB,muZ] = FEmat_dipoles(mesh,massCoeff,q)
% Constructs the dipole-moment matrices for spatially-invariant electric or magnetic fields.
% The output will be a cell array, where each entry corresponds to the dipole
% operator for the field along each axis.
arguments
	mesh FE_mesh
	massCoeff {mustBeA(massCoeff,{'FE_coeff','function_handle','char','string','numeric'})} = [];
	q (1,1) double {mustBeNumeric,mustBeFinite} = -1;
end
[numDim] = mesh.num("dim");

% compute electric dipole moments
muE = cell(1,numDim);
for ii = 1:numDim
	E = FE_coeff(mesh,-mesh.ptsPoly(ii,:));
	muE{ii} = q * FEmat(mesh,"scalar",E);
end

% exit without computing 1st/2nd-order magnetic dipole moments
if nargout <= 1; return; end

switch numDim
	case 2
		% vector potential for:
		Az = FE_coeff(mesh,"@(x,y) -1/2 * [y; -x]"); % z-aligned B-field

		% 1st-order magnetic dipole moment
		muB{1} = 1i*q * FEmat(mesh,"skew-vector",massCoeff,Az);

		% exit without computing 2nd-order magnetic dipole moment
		if nargout <= 2; return; end

		% 2nd-order magnetic dipole moment
		muZ{1} = q^2 * FEmat(mesh,"symm-vector",massCoeff,Az,Az);

	case 3
		% vector potentials for:
		A{1} = FE_coeff(mesh,"@(x,y,z) -1/2 * [0*x; z; -y]"); % x-aligned B-field
		A{2} = FE_coeff(mesh,"@(x,y,z) -1/2 * [-z; 0*y; x]"); % y-aligned B-field
		A{3} = FE_coeff(mesh,"@(x,y,z) -1/2 * [y; -x; 0*z]"); % z-aligned B-field
		
		% 1st-order magnetic dipole moments
		muB = cell(1,numDim);
		for ii = 1:numDim
			muB{ii} = 1i*q * FEmat(mesh,"skew-vector",massCoeff,A{ii});
		end

		% exit without computing 2nd-order magnetic dipole moments
		if nargout <= 2; return; end

		% 2nd-order magnetic dipole moments
		muZ = cell(numDim,numDim);
		for ii = 1:numDim
			for jj = 1:numDim
				muZ{ii,jj} = q^2 * FEmat(mesh,"symm-vector",massCoeff,A{ii},A{jj});
			end
		end
end
% DONE
end