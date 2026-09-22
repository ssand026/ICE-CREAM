function [sphereC,sphereR,newPts,newTri,newReg] = enclosing_sphere(pts,tri,reg,opt)
% generates a spherical region around some initial mesh object
arguments
	pts (:,:) double {mustBeReal,mustBeFinite}
	tri (:,:) double {mustBeInteger,mustBePositive}
	reg (1,:) double {mustBeInteger,mustBeNonnegative} = [];
	opt.radiusScale (1,1) double {mustBeGreaterThanOrEqual(opt.radiusScale,1)} = 1.02;
	opt.exteriorGap (1,1) double {mustBeGreaterThanOrEqual(opt.exteriorGap,1)} = 1.5;
	opt.smoothMesh (1,1) logical = true;
end

% algorithm parameters
numDim = size(pts,1);
MAX_ITERS = 20;
EDGE_SIGMA = 1.0;
THRESH_QUAL = 1/3;
THRESH_VOL = numDim+1;
useConstrainedDT = true; % my constrained delaunay-triangulation free approach is not fully working 

%% get mesh properties
%===========================================================
if (numDim ~= 2) && (numDim ~= 3)
	error("ERROR: the number of spatial dimensions must be either 2 or 3")
end

% get the set of nodes/faces that form the domain boundary
warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
intTR = triangulation(tri.', pts.');
bndFaces = freeBoundary(intTR);
bndNodes = unique(bndFaces(:));
bndPtsMap = accumarray(bndNodes(:),(1:numel(bndNodes)).',[max(bndNodes),1]);
bndFaces = sortrows(sort(bndPtsMap(bndFaces),2));

% get the edge lengths of all boundary faces
bndPts = intTR.Points(bndNodes,:);
bndEdges = subfaces(bndFaces.',2);
bndEdgeLen = hypervolumes(bndPts.',bndEdges.');
numBndPts = numel(bndNodes);

%% find the enclosing sphere center and radius
%===========================================================
% find the smallest sphere that encloses the domain
[outerC,outerR] = outer_sphere(bndPts.');

% set the target edge length for the exterior surface
targLen = mean(bndEdgeLen) + EDGE_SIGMA * std(bndEdgeLen);

% get the volume of a regular simplex with side length equal to targLen
simplexVol = sqrt((numDim+1)/(2^numDim))/factorial(numDim) * (targLen^numDim);

% determine the enclosing sphere radius
simplexHeight = sqrt((numDim+1)/(2*numDim)); % height of a simplex with an edge length of 1
sphereR = max(opt.radiusScale * outerR, outerR + opt.exteriorGap * (1+simplexHeight)/2 * targLen);
sphereC = outerC;

% exit if mesh is not needed
if nargout <= 2
	newPts = [];
	newTri = [];
	newReg = [];
	return
end

%% generate points on the enclosing sphere
%===========================================================
sphereSurfVol = @(r,dim) 2*pi/factorial(dim-1) * (r^(dim-1));
numSpherePts = ceil(sphereSurfVol(sphereR,numDim)/(targLen^(numDim-1)));

% ensure that the number of boundary points is odd
numSpherePts = numSpherePts + mod(numSpherePts,2) + 1;

switch numDim
	case 2
		spherePts = linspace(0,2,numSpherePts+1);
		spherePts = [cospi(spherePts(1:end-1)); sinpi(spherePts(1:end-1))].';
		spherePts = (sphereR * spherePts) + sphereC;
	case 3
		spherePts = (sphereR * nsphere(numSpherePts)) + sphereC;
end


%% refine the exterior-region mesh
%===========================================================
% Working on a method that does not require the external 3-D constrained
% Delaunay triangulation algorithm, leaving it here despite not working
if useConstrainedDT
	%% combine the domain boundary and exterior boundary points
	%===========================================================
	% find the largest sphere that is enclosed by the domain
	[innerC,~] = inner_sphere(bndPts.',intTR);
	
	% invert boundary coordinates around the inscribed midpoint
	invPts = invert_pts([bndPts; spherePts],innerC);
	invPts(end+1,:) = 0; % append point at "infinity"
	
	% perform delaunay triangulation on the inverted coordinates. The domain
	% boundary now lies on the mesh exterior, allowing the faces to be constrained.
	invTR = constrained_DT(invPts,bndFaces);
	
	% remove triangulations containing any non-boundary points
	invTri = invTR.ConnectivityList;
	invTri(all(isbetween(invTri-numBndPts,1,numSpherePts+1),2),:) = [];
	invTR = triangulation(invTri, invPts);
	
	% un-invert the inverse triangulation to obtain the triangulated exterior points
	extPts = revert_pts(invTR.Points,innerC);
	extTri = invTR.ConnectivityList;
	infPts = find(any(~isfinite(extPts),2));
	infTri = any(ismember(extTri,infPts),2);
	extPts(infPts,:) = 0;
	extTri(infTri,:) = [];
	extTR = triangulation(extTri,extPts);
	
	%% add points within the exterior region
	%===========================================================
	for iter = 1:MAX_ITERS
		% find triangles in the exterior mesh that meet the criteria to encompass a new point
		triQual = triqual(extTR.Points.',extTR.ConnectivityList.');
		triVol = hypervolumes(extTR.Points.',extTR.ConnectivityList.');
		% inInterior = ~isnan(pointLocation(intTR,incenter(extTR)));
		inInterior = false; % ^^^ should always return false in this setup
		validTri = find((~inInterior) & (triVol(:) > THRESH_VOL*simplexVol));
		idealTri = validTri(triQual(validTri) > THRESH_QUAL);
		
		if isempty(validTri)
			% the exterior region has been fully meshed
			break
		elseif isempty(idealTri)
			% add a single point to the exterior region
			[~,bestTri] = max(triQual(validTri));
			newPts = tricenters(extTR,validTri(bestTri));
		else
			% add new points to the exterior region
			newPts = tricenters(extTR,idealTri);
		end
		
		% add the new points to the inverted mesh
		newPts = invert_pts(newPts,innerC);
		invTR = constrained_DT([invTR.Points; newPts],bndFaces);
		
		% remove outside points
		invTri = invTR.ConnectivityList;
		invTri(all(isbetween(invTri-numBndPts,1,numSpherePts+1),2),:) = [];
		invTR = triangulation(invTri, invTR.Points);
		
		% revert to non-inverted coordinates
		extPts = revert_pts(invTR.Points,innerC);
		extTri = invTR.ConnectivityList;
		infPts = find(any(~isfinite(extPts),2));
		infTri = any(ismember(extTri,infPts),2);
		extPts(infPts,:) = 0;
		extTri(infTri,:) = [];
		extTR = triangulation(extTri,extPts);
	end
	
	% remove interior/infinite triangulations
	extPts = extTR.Points;
	extTri = extTR.ConnectivityList;
	extPtsSubset = setdiff(1:size(extPts,1),infPts);
	extTriSubset = all(ismember(extTri,extPtsSubset),2);
	extPtsMap = accumarray(extPtsSubset(:),(1:numel(extPtsSubset)).',[max(extPtsSubset),1]);
	extPts = extPts(extPtsSubset,:);
	extTri = extPtsMap(extTri(extTriSubset,:));
	extTR = triangulation(extTri,extPts);
	
	if numDim == 2
		% remove triangles in extTR that are enclosed within the region intTR
		inInterior = ~isnan(pointLocation(intTR,incenter(extTR)));
		extTri = extTR.ConnectivityList;
		extTri(inInterior,:) = [];
		extTR = triangulation(extTri,extTR.Points);
	end
else
	%% combine the points on the domain boundary and sphere
	%===========================================================
	% perform initial triangulation
	extDT = delaunayTriangulation([bndPts; spherePts]);
	
	% remove triangulation of the interior region
	inOrig = ~isnan(pointLocation(intTR,tricenters(extDT)));
	extTri = extDT.ConnectivityList;
	extTri(inOrig,:) = [];
	extTR = triangulation(extTri,extDT.Points);
	
	%% add points within the exterior region
	%===========================================================
	for iter = 1:MAX_ITERS
		% find triangles in the exterior mesh that meet the criteria to encompass a new point
		triQual = triqual(extTR.Points.',extTR.ConnectivityList.');
		triVol = hypervolumes(extTR.Points.',extTR.ConnectivityList.');
		inInterior = ~isnan(pointLocation(intTR,incenter(extTR)));
		validTri = find((~inInterior) & (triVol(:) > THRESH_VOL*simplexVol));
		idealTri = validTri(triQual(validTri) > THRESH_QUAL);
		
		if isempty(validTri)
			% the exterior region has been fully meshed
			break
		elseif isempty(idealTri)
			% add a single point to the exterior region
			[~,bestTri] = max(triQual(validTri));
			newPts = tricenters(extTR,validTri(bestTri));
		else
			% add new points to the exterior region
			newPts = tricenters(extTR,idealTri);
		end
		
		% add the new points to the mesh
		extDT = delaunayTriangulation([extTR.Points; newPts]);
		inOrig = ~isnan(pointLocation(intTR,tricenters(extDT)));
		
		% remove triangulation of the interior region
		extTri = extDT.ConnectivityList;
		extTri(inOrig,:) = [];
		extTR = triangulation(extTri,extDT.Points);
	end
	
	%% resolve conflicts between the internal/external regions
	%===========================================================
	for iter = 1:MAX_ITERS
		% extract portion of the triangulation next to a boundary face
		extTri = extTR.ConnectivityList;
		triNearAbsent = any(extTri<=numBndPts,2);
		nearbyTri = extTri(triNearAbsent,:);
		
		% check if any boundary faces of the internal region are missing from the external mesh
		extBnds = subfaces(nearbyTri.',numDim);
		extBnds = extBnds(all(extBnds<=numBndPts,2),:);
		
		isAbsent = ~ismember(bndFaces,extBnds,"rows");
		absentFaces = bndFaces(isAbsent,:);
		
		if nnz(isAbsent) <= 0
			break
		end
		
		% find shared edges in the missing face list
		edgePerms = ncombsk(numDim,2).';
		absentEdges = reshape(absentFaces(:,edgePerms(:)).',2,[]).';
		sharedEdges = unique(absentEdges(~isunique(absentEdges,1),:),"rows");
		
		% restore missing boundary faces by searching over the shared edges
		addTri = [];
		isFixed = false(size(absentFaces,1),1);
		toRemove = false(size(nearbyTri,1),1);
		for ii = 1:size(sharedEdges,1)
			% find triangles in contact with the faces containing the shared edge
			absEdge = sharedEdges(ii,:);
			edgeInFace = (sum(ismember(absentFaces,absEdge),2) >= 2);
			absFace = absentFaces(edgeInFace,:);
			faceInTri = (sum(ismember(nearbyTri,absFace(:)),2) >= numDim);
			adjTri = nearbyTri(faceInTri,:);
			
			% find node/s in the adjacent tri that do not lie on the faces
			other = setdiff(adjTri,absFace);
			other = other(~ismember(other,bndNodes));
			
			% perform fixes depending on configuration
			if nnz(faceInTri) < 2
				% there are too few triangles
				% adjTri = nearbyTri((sum(ismember(nearbyTri,absFace(:)),2) >= 2),:)
			elseif nnz(faceInTri) > 2
				% there are too many triangles
			elseif isempty(other)
				% add a new simplex to restore missing faces
				addTri(end+1,:) = unique(absFace(:)).';
				isFixed(edgeInFace) = true;
			elseif isscalar(other)
				% swap quad edges to restore missing faces
				nearbyTri(faceInTri,:) = [absFace,[other; other]];
				isFixed(edgeInFace) = true;
			elseif numel(other)==2
				% try to perform swap on whole section
				oppEdge = intersect(adjTri(1,:),adjTri(2,:));
				isThird = all(ismember(nearbyTri,[oppEdge(:); other(:)]),2);
				diagPnt = adjTri(any(adjTri==other(2),2) & ismember(adjTri,absEdge));
				if nnz(isThird)==1 && numel(oppEdge)==2
					addTri(end+1:end+4,:) = [...
						other(1),absFace(1,:);
						other(1),absFace(2,:);
						other(1),other(2),diagPnt,oppEdge(1);
						other(1),other(2),diagPnt,oppEdge(2)];
					% nearbyTri(isThird | edgeInTri,:) = NaN;
					toRemove = toRemove | isThird | faceInTri;
					isFixed(edgeInFace) = true;
				end
			end
		end
		
		% take care of the remaining faces
		absentFaces = absentFaces(~isFixed,:);
		for ii = 1:size(absentFaces,1)
			% get edges of the triangulation where one node lies on the missing face
			absFace = absentFaces(ii,:);
			adjTri = nearbyTri((sum(ismember(nearbyTri,absFace),2) >= numDim-1),:);
			adjEdg = subfaces(adjTri.',2);
			onFace = ismember(adjEdg,absFace);
			
			isAdjacent = xor(onFace(:,1),onFace(:,2));
			adjEdg = adjEdg(isAdjacent,:);
			onFace = onFace(isAdjacent,:);
			
			% find the node common to all connected edges
			common = find(sum(sparse(adjEdg(onFace),adjEdg(~onFace),1),1)>=numDim);
			common = common((common > numBndPts));
			if isscalar(common)
				% add the missing simplex to the triangulation
				addTri(end+1,:) = [absFace,common];
			else
				%???
			end
		end
		
		% update the triangulation
		extTri(triNearAbsent,:) = nearbyTri;
		% extTri(any(isnan(extTri),2),:) = [];
		toRmv = find(triNearAbsent);
		toRmv = toRmv(toRemove);
		extTri(toRmv,:) = [];
		extTR = triangulation([extTri;addTri],extTR.Points);
	end
	
end

%% combine the interior and exterior meshes
%===========================================================
numIntPts = size(intTR.Points,1);
numExtPts = size(extTR.Points,1);

numIntTri = size(intTR.ConnectivityList,1);
numExtTri = size(extTR.ConnectivityList,1);

% map the original boundary nodes indices to their new locations
extPtsOld = [(1:numExtPts).'];
extPtsNew = [bndNodes(:); numIntPts+(1:(numExtPts-numBndPts)).'];
extPtsMap = full(sparse(extPtsOld,1,extPtsNew,max(extPtsOld),1,max(extPtsNew)));

% append exterior mesh to the interior mesh
newPts = [intTR.Points; extTR.Points(numBndPts+1:end,:)].';
newTri = [intTR.ConnectivityList; extPtsMap(extTR.ConnectivityList)].';
if isempty(reg)
	newReg = [zeros(1,numIntTri),ones(1,numExtTri)];
elseif isscalar(reg)
	newReg = [(reg)*ones(1,numIntTri),(reg+1)*ones(1,numExtTri)];
else
	newReg = [reg, (max(reg)+1)*ones(1,numExtTri)];
end

% apply smoothing
if opt.smoothMesh
	newPts = mesh_smooth(newPts,newTri,newReg,"fixedSurfs","all","iters",20);
	newPts(:,1:numIntPts) = intTR.Points.';
end
% DONE
end

% external functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TR_mid] = tricenters(TR,index)
% get the simplex midpoints for a triangulation
arguments
	TR
	index (:,1) = []
end

if isempty(index)
	tri = TR.ConnectivityList;
else
	tri = TR.ConnectivityList(index,:);
end

numDim = size(TR.Points,2);
[numTri,numVtx] = size(tri);

TR_mid = zeros(numTri,numDim);
for ii = 1:numVtx
	vtx_ii = tri(:,ii);
	TR_mid = TR_mid + TR.Points(vtx_ii,:);
end
TR_mid = TR_mid/numVtx;
% DONE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DT] = constrained_DT(points,constraint)
% performs a constrained delaunay triangulations
arguments
	points (:,:) double
	constraint (:,:) double
end
numDim = size(points,2);
if numDim == 2
	% use the builtin matlab method
	DT = delaunayTriangulation(points,constraint);
elseif numDim == 3
	% try calling the GibbonCode method
	try
		[~,DT] = evalc("constrainedDelaunayTetGen(points,constraint)");
	catch ME
		if strcmp(ME.identifier,"MATLAB:UndefinedFunction")
			msg = ...
				"The function 'constrainedDelaunayTetGen' is not found on the current path. " + ...
				"Try installing the GibbonCode package for MATLAB to use its implementation " + ...
				"of a constrained 3D Delaunay-triangulation algorithm:" + newline + ME.message;
			throw(MException(ME.identifier,msg))
		else
			throw(ME)
		end
	end
else
	error("ERROR: the point set for contrained_DT must be two or three-dimensional")
end
% DONE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [invPts] = invert_pts(origPts,center)
% invert point locations across the specified center point
ptsDist = sum((origPts-center).^2,2);
invPts = (origPts-center) ./ ptsDist;
% DONE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [origPts] = revert_pts(invPts,center)
% un-invert point locations across the specified center point
ptsDist = sum((invPts).^2,2);
origPts = (invPts ./ ptsDist) + center;
% DONE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pts] = nsphere(numPts)
% creates a psuedo-evenly spaced set of points on the unit sphere using the
% Fibonacci spiral method
arguments
	numPts (1,1) double {mustBePositive,mustBeInteger}
end
theta = pi*(sqrt(5)-1).*(0:numPts-1);
z = linspace(-1,+1,numPts);
r = sqrt(1 - z.^2);
x = r .* cos(theta);
y = r .* sin(theta);
pts = [x(:),y(:),z(:)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%