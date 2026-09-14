function [muB,muZ] = magnetic_dipole(mesh,mass,q)
% Constructs the dipole-moment matrices for spatially-invariant magnetic fields.
% The first output will be a cell array, where each entry corresponds to the 
% first-order magnetic dipole operator for the field along each axis.
% The second output will be a cell array, where each entry corresponds to the
% second-order magnetic dipole operator for combination of fields along each axis.
arguments
	mesh FE_mesh
	mass {mustBeA(mass,{'FE_coeff','function_handle','char','string','numeric'})} = 1;
	q (1,1) double {mustBeNumeric,mustBeFinite} = -1;
end

[numDim] = mesh.num("dim");

switch numDim
	case 2
		% vector potential for:
		Az = FE_coeff(mesh,"@(x,y) -1/2 * [y; -x]"); % z-aligned B-field

		% 1st-order magnetic dipole moment
		muB{1} = 1i*q * FEmat(mesh,"skew-vector",mass,Az);

		% exit without computing 2nd-order magnetic dipole moment
		if nargout <= 1; return; end

		% 2nd-order magnetic dipole moment
		muZ{1} = q^2 * FEmat(mesh,"symm-vector",mass,Az,Az);

	case 3
		% vector potentials for:
		A{1} = FE_coeff(mesh,"@(x,y,z) -1/2 * [0*x; z; -y]"); % x-aligned B-field
		A{2} = FE_coeff(mesh,"@(x,y,z) -1/2 * [-z; 0*y; x]"); % y-aligned B-field
		A{3} = FE_coeff(mesh,"@(x,y,z) -1/2 * [y; -x; 0*z]"); % z-aligned B-field
		
		% 1st-order magnetic dipole moments
		muB = cell(1,numDim);
		for ii = 1:numDim
			muB{ii} = 1i*q * FEmat(mesh,"skew-vector",mass,A{ii});
		end

		% exit without computing 2nd-order magnetic dipole moments
		if nargout <= 1; return; end

		% 2nd-order magnetic dipole moments
		muZ = cell(numDim,numDim);
		for ii = 1:numDim
			for jj = 1:numDim
				muZ{ii,jj} = q^2 * FEmat(mesh,"symm-vector",mass,A{ii},A{jj});
			end
		end
end
% DONE
end