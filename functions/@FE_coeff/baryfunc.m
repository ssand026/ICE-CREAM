function [bf] = baryfunc(mesh)
% generates a FE_coeff object corresponding to the generic
% set of barycentric basis functions for the given mesh
arguments
	mesh FE_mesh
end
[numTri,numPoly] = mesh.num("tri","poly");
coeffs = ones(1,1,numPoly,numTri);
bf = FE_coeff(mesh,coeffs,mesh.vtxPoly,(1:numPoly).');
% DONE
end