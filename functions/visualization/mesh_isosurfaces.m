function [surfNodes, surfFaces] = mesh_isosurfaces(pts,tri,data,vals,opt)
% MESH_ISOSURFACES Given a scalar field defined on some mesh, returns the set of nodes/faces that
% correspond to the isosurfaces with the requested values.
% 
% SEE ALSO: MESH_ISOSURFACE
arguments
	pts (:,:) double {mustBeReal,mustBeFinite}
	tri (:,:) {mustBeInteger,mustBePositive}
	data (:,:) double {mustBeReal,mustBeFinite}
	vals (1,:) double {mustBeReal,mustBeFinite}
	opt.tolerance (1,1) {mustBeBetween(opt.tolerance,0,1)} = 1e-9;
end

numSurf = numel(vals);
surfNodes = cell(1,numSurf);
surfFaces = cell(1,numSurf);

% check data size
%-----------------------------------------------------------
[~,numPts] = size(pts);
[~,numTri] = size(tri);
if isscalar(data)
	% data is constant and continuous
	return
elseif isvector(data) && numel(data)==numTri
	% data is constant and disjoint
	return
elseif isvector(data) && numel(data)==numPts
	% data is variable and continuous
	data = data(tri);
elseif isequal(size(data),size(tri))
	% data is variable and disjoint
else
	error("ERROR: data has invalid shape")
end

% ignore regions of the mesh where the data is invariant
%-----------------------------------------------------------
tol = max(opt.tolerance/eps,1) * max(eps(data),[],1);
isFlat = ((max(data,[],1) - min(data,[],1)) <= tol);
tri = tri(:,~isFlat);
data = data(:,~isFlat);
tol = tol(:,~isFlat);

% ignore regions of the mesh where the field does not span any isosurface
%-----------------------------------------------------------
sortedVals = sort(vals);
dataMax = max(data,[],1) - tol;
dataMin = min(data,[],1) + tol;
noSpan = (dataMax < min(vals)) | (dataMin > max(vals));
if all(noSpan)
	% data does not span any surface value
	return
end
for ii = 2:numSurf
	above = (dataMin > sortedVals(ii-1));
	below = (dataMax < sortedVals(ii));
	noSpan = noSpan | (above & below);
	if all(noSpan)
		% data does not span any surface value
		return 
	end
end
tri = tri(:,~noSpan);
data = data(:,~noSpan);

% find the nodes/faces that form each isosurface
%-----------------------------------------------------------
for ii = 1:numSurf
	[surfNodes{ii}, surfFaces{ii}] = mesh_isosurface(pts,tri,data,vals(ii));
end

% DONE
end