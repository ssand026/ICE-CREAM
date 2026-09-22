function [nodes, faces] = mesh_isosurface(pts,tri,data,value)
% MESH_ISOSURFACE Given a scalar field defined on some mesh, returns the set of nodes/faces that
% correspond to thre requested isosurface value.
arguments
	pts (:,:) double {mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	data (:,:) double {mustBeReal,mustBeFinite}
	value (1,1) double {mustBeFinite} = 0;
end

% circumvent the costly overhead of "unique(...)" using:
%	[C,ia,ic] = matlab.internal.math.uniquehelper(A,doSort,isFirst,byRows)
fast_unique = @(x,opts) matlab.internal.math.uniquehelper(x,opts(1),opts(2),opts(3));

% subtract the surface-value from the data, we now just need
% to find the zero-valued surface
if (value ~= 0)
	data = data - value;
end

%% check size of data
%===========================================================
[numDim,numPts] = size(pts);
[numVtx,numTri] = size(tri);
if isscalar(data) || isempty(data)
	% data is constant and continuous
	nodes = [];
	faces = [];
	return
elseif isvector(data) && (numel(data)==numTri)
	% data is constant and disjoint
	nodes = [];
	faces = [];
	return
elseif isvector(data) && (numel(data)==numPts)
	% data is variable and continuous
	data = data(tri);
elseif isequal(size(data),size(tri))
	% data is variable and disjoint
else
	error("ERROR: data has invalid shape")
end

%% find any existing faces that lie on the zero-valued surface
%===========================================================
% check if the data contains any zero-valued points
[dataMin,dataMax] = bounds(data(:));
if (dataMin > 0) || (dataMax < 0)
	% data does not intersect with the zero-valued surface
	nodes = [];
	faces = [];
	return
end

% determine if the simplex contains ...
%-----------------------------------------------------------
exact = logical(data == 0);
numExact = sum(exact,1);
oneZero = (numExact == (numVtx-1)); % one zero-valued face
allZero = (numExact >= (numVtx+0)); % all zero-valued faces

% extract any existing zero-valued nodes/faces
%-----------------------------------------------------------
zeroSurf = [];
hasZero = (oneZero | allZero);
if any(hasZero)
	% extract the existing zero-valued faces
	if any(oneZero)
		zeroSurf = [zeroSurf, reshape(tri(exact & oneZero),numVtx-1,[])];
	end
	if any(allZero)
		facePerms = ncombsk(numVtx,numVtx-1).';
		zeroSurf = [zeroSurf, reshape(tri(facePerms(:),allZero),numVtx-1,[])];
	end

	% combine zero-face data with node positions
	zeroSurf = reshape(pts(:,zeroSurf(:)),numDim,numVtx-1,[]);
	zeroSurf = permute(zeroSurf,[2,3,1]);

	% ignore zero-faced simplexes going forward
	tri = tri(:,~hasZero);
	data = data(:,~hasZero);
end

%% interpolate over each simplex to find the zero-valued surfaces
%===========================================================
% check for zero-crossings
if (dataMin >= 0) || (dataMax <= 0)
	% data does not span the zero-valued surface
	doInterpolation = false;
else
	% check if each simplex contains a zero-crossing
	above = logical(data > 0);
	below = logical(data < 0);
	spansZero = any(above,1) & any(below,1);

	% ignore non-spanning simplexes going forward
	above = above(:,spansZero);
	below = below(:,spansZero);
	data = data(:,spansZero);
	tri = tri(:,spansZero);
	[~,numTri] = size(tri);

	% determine if any data spans the zero-valued surface
	doInterpolation = any(spansZero(:));
end

if doInterpolation
	% interpolate data to find the zero-valued surface
	%-----------------------------------------------------------
	% get the set of all simplex edges
	[edgeA, edgeB] = find(triu(ones(numVtx)));
	numEdge = numel(edgeA & edgeB);
	
	% find the set of orthogonal edge pairs
	orthEdges = (edgeA'~=edgeB) & (edgeB'~=edgeA) & (edgeA'~=edgeA) & (edgeB'~=edgeB);
	[row,col] = find(orthEdges);
	orthPairs = fast_unique(sort([row,col],2),[false,true,true]);
	
	% for each set of simplex edges, find the zero point location, if it exists
	edges = false(numEdge,numTri);
	zeroPts = NaN(numDim,numTri,numEdge);
	
	for ii = 1:numEdge
		isAB = above(edgeA(ii),:) & below(edgeB(ii),:);
		isBA = above(edgeB(ii),:) & below(edgeA(ii),:);
		edges(ii,:) = xor(isAB,isBA);
		
		edgeFrac = data([edgeB(ii),edgeA(ii)],edges(ii,:));
		edgeFrac = abs(edgeFrac) ./ abs(edgeFrac(1,:) - edgeFrac(2,:));
	
		ptsA = tri(edgeA(ii),edges(ii,:));
		ptsB = tri(edgeB(ii),edges(ii,:));
		midpt = (edgeFrac(1,:) .* pts(:,ptsA)) + (edgeFrac(2,:) .* pts(:,ptsB));
	
		zeroPts(:,edges(ii,:),ii) = midpt;
	end
	zeroPts = permute(zeroPts,[3,2,1]);
	numCross = sum(edges,1);
	
	% get the subset of zeroPts where ...
	%-----------------------------------------------------------
	% the interpolated zero-surface is a simplex
	setA = edges(:,(numCross==(numVtx-1)));
	zeroPtsA = zeroPts(:,(numCross==numVtx-1),:);
	surfA = zeroPtsA(repmat(setA,1,1,numVtx-1));
	surfA = reshape(surfA,numVtx-1,[],numVtx-1);
	
	% the interpolated zero-surface is NOT a simplex
	setB = edges(:,(numCross==numVtx));
	setC = edges(:,(numCross==numVtx));
	zeroPtsBC = zeroPts(:,(numCross==numVtx),:);
	
	whereOrth = setB(orthPairs(:,1),:) & setC(orthPairs(:,2),:);
	whereOrth = max((1:size(orthPairs,1)).' .* whereOrth,[],1);
	
	b = orthPairs(whereOrth,1).';
	setB(sub2ind(size(setB),b,1:numel(b))) = false;
	surfB = zeroPtsBC(repmat(setB,1,1,numVtx-1));
	surfB = reshape(surfB,numVtx-1,[],numVtx-1);
	
	c = orthPairs(whereOrth,2).';
	setC(sub2ind(size(setC),c,1:numel(c))) = false;
	surfC = zeroPtsBC(repmat(setC,1,1,numVtx-1));
	surfC = reshape(surfC,numVtx-1,[],numVtx-1);

	% combine all zero-surfaces
	zeroSurf = [zeroSurf, surfA, surfB, surfC];
end

%% extract the unique nodes/faces from the surface array
%===========================================================
if isempty(zeroSurf)
	nodes = [];
	faces = [];
else
	zeroSurf = reshape(zeroSurf,[],numVtx-1);
	[nodes,~,faces] = fast_unique(zeroSurf,[false,true,true]);
	faces = reshape(faces,numVtx-1,[]).';
	faces = fast_unique(sort(faces,2),[false,true,true]); %#ok
end

%% DONE
end