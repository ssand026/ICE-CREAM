function [intoB,exitB] = femat_dirichlet(mesh,nodes)
% Finite-element basis-reduction/expansion matrices for Dirichlet boundary-conditions
%
% The output 'intoB' is the basis reduction matrix
% The output 'exitB' is the basis expansion matrix
%
% Calculations done in the reduced basis will automatically satisfy Dirichlet
% boundary conditions on the specified nodes.
arguments
	mesh FE_mesh
	nodes (1,:) {mustBeInteger,mustBePositive} = mesh.boundaryNodes;
end
numPts = mesh.num("pts");
onBnd = (sum(mesh.vtxConn(nodes,:),1) == sum(mesh.vtxConn,1));

% get the dirichlet boundary-condition application matrix
intoB = speye(numPts,numPts);
intoB(onBnd,:) = [];

% get the dirichlet boundary-condition reversal matrix
exitB = intoB.';

% DONE
end