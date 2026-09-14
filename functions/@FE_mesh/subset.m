function [newMesh] = subset(mesh,param,index)
% Returns a mesh corresponding to a subset of the original mesh
arguments
	mesh FE_mesh
	param (1,1) {mustBeMember(param,["pts","tri","reg"])}
	index (:,1) {mustBeNumericOrLogical}
end

% parse the sub-section indices
if isnumeric(index)
	if any(index < 0) || (~isinteger(index) && any(mod(index,1)~=0,"all"))
		% indices are not positive integers
		error("ERROR: the elements of 'index' must be non-negative integers")
	elseif any(index > mesh.num(param))
		% indices exceeds array bounds
		error("ERROR: elements of 'index' exceed the number of "+param)
	elseif any(~isunique(index))
		% index contains repeats
		error("ERROR: the input 'index' contains repeated entries")
	elseif param == "reg" && ~isempty(setdiff(index,unique(mesh.reg)))
		error("ERROR: the input 'index' references non-existent regions")
	end
elseif islogical(index)
	if numel(index) ~= mesh.num(param)
		% index exceeds array bounds
		error("ERROR: length of 'index' does not equal the number of "+param)
	else
		% convert to a numeric index
		index = find(index);
	end
end

% extract portion of original mesh
switch param
	case "pts"
		ptsSubset = index;
		triSubset = all(ismember(mesh.tri,index),1);
	case "tri"
		triSubset = index;
		ptsSubset = unique(mesh.tri(:,triSubset));
	case "reg"
		triSubset = ismember(mesh.reg,index);
		ptsSubset = unique(mesh.tri(:,triSubset));
end

% remove unused pts/tri and re-assign pts labels
ptsMap = accumarray(ptsSubset,(1:numel(ptsSubset)).',[max(ptsSubset),1]);
newPts = mesh.pts(:,ptsSubset);
newTri = ptsMap(mesh.tri(:,triSubset));

% take region subset
if isscalar(mesh.reg)
	newReg = mesh.reg;
else
	newReg = mesh.reg(:,triSubset);
end

% construct a new mesh
newMesh = FE_mesh(newPts,newTri,newReg,"degree",mesh.degree);
% DONE
end