function [obj] = restructure(obj,method)
% Re-orders the mesh to minimize the fill-in of the decomposed adjacency matrix.
% Will generally improve the performance of calculations with large/sparse systems.
arguments
	obj FE_mesh
	method (1,1) {mustBeMember(method,["amd","symamd","symrcm","dissect"])}
end

% get the current sparse-solver parameters
spvals = spparms;

% configure spparams for better ordering outcomes
spparms("tight");

% select ordering method
switch method
	case "amd"
		% use the Approximate Minimum Degree ordering (best for iterative solvers)
		ordering = amd(obj.adjacencyMatrix);

	case "symamd"
		% use the symmetric Approximate Mininimum Degree ordering
		ordering = symamd(obj.adjacencyMatrix);

	case "symrcm"
		% use the symmetric Reverse Cuthill-McKee ordering algorithm
		ordering = symrcm(obj.adjacencyMatrix);
	
	case "dissect"
		% use the dissection ordering algorithm (best for the default '\' solver)
		ordering = dissect(obj.adjacencyMatrix,"NumIterations",50);
end

% re-order mesh
obj = obj.reorder("pts",ordering);

% restore the previous sparse-solver parameters
spparms(spvals);
% DONE
end