function [obj] = reorder(obj,param,ordering)
% Re-structures the mesh according to the specified node ordering
arguments
	obj FE_mesh
	param (1,1) {mustBeMember(param,["pts","tri"])}
	ordering (1,:) {mustBeInteger,mustBePositive}
end

% list where each element contains the corresponding node's new location
[sortedOrdering,orderingMap] = sort(ordering);

switch param
	case "pts"
		% reorder the pts array
		if ~isequal(sortedOrdering,[1:obj.num("pts")])
			error("ERROR: the input 'ordering' is not a permutation of the node indices")
		end
		
		% update the mesh
		obj.pts = obj.pts(:,ordering);
		obj.tri = orderingMap(obj.tri);

	case "tri"
		% reorder the tri/reg array
		if ~isequal(sortedOrdering,[1:obj.num("tri")])
			error("ERROR: the input 'ordering' is not a permutation of the triangulation indices")
		end

		% update the mesh
		obj.tri = obj.tri(:,ordering);
		if ~isscalar(obj.reg)
			obj.reg = obj.reg(:,ordering);
		end
end
% DONE
end