function [newPts,newTri] = mesh_subset(pts,tri,param,index)
% returns a new sub-mesh corresponding to 
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	param (1,1) {mustBeMember(param,["pts","tri"])}
	index (:,1) {mustBeA(index,{'numeric','logical'})}
end
numPts = size(pts,2);
numTri = size(tri,2);

% parse the sub-section indices
if isnumeric(index)
	if any(index < 1) || (~isinteger(index) && any(mod(index,1)~=0,"all"))
		% index elements are not positive integers
		error("ERROR: the elements of 'index' must be positive integers")
	elseif ((param=="pts") && any(index > numPts)) || ((param=="tri") && any(index > numTri))
		% index elements exceed the array bounds
		error("ERROR: elements in 'index' exceed the number of " + param)
	elseif any(~isunique(index))
		% index contains repeats
		error("ERROR: the input 'index' contains repeated entries")
	end
elseif islogical(index)
	if ((param=="pts") && (numel(index) ~= numPts)) || ((param=="tri") && (numel(index) ~= numTri))
		% index length exceeds array bounds
		error("ERROR: length of 'index' does not equal the number of "+param)
	else
		% convert to a numeric index
		index = find(index);
	end
end

% find indices that extract the relevant portion of original mesh
index = cast(index,"like",tri);
switch param
	case "pts"
		ptsSubset = index;
		triSubset = all(ismember(tri,index),1);
	case "tri"
		ptsSubset = unique(tri(:,index));
		triSubset = index;
end

% remove unused pts/tri and re-assign node labels
ptsMap = accumarray(ptsSubset,(1:numel(ptsSubset)).',[max(ptsSubset),1]);
newPts = pts(:,ptsSubset);
newTri = ptsMap(tri(:,triSubset));

% DONE
end