function [obj] = load(geom,nodes,opt)
% Generates a FE_mesh from the given input.
% The resulting mesh will approximately contain the requested number of nodes.
%
% REQUIRES:
%	Partial Differential Equation Toolbox
%
% SEE ALSO: FEGEOMETRY, GENERATEMESH, DECSG
arguments
	geom {mustBeA(geom,["string","char","function_handle","numeric","FEMesh","fegeometry", ...
		"triangulation","delaunayTriangulation","pde.DiscreteGeometry","pde.AnalyticGeometry"])}
	nodes (1,1) double {mustBeInteger,mustBePositive}
	opt.tolerance (1,1) double {mustBeInRange(opt.tolerance,0,1)} = 1e-3;
	opt.converged (1,1) double {mustBeInRange(opt.converged,0,1)} = 1e-5;
	
	opt.Hgrad (1,1) double {mustBeInRange(opt.Hgrad,1,2)} = sqrt(2);
	opt.maxIters double {mustBeScalarOrEmpty,mustBeInteger,mustBeNonnegative} = [];
	opt.normalize logical {mustBeScalarOrEmpty} = [];
	opt.cullEdges logical {mustBeScalarOrEmpty} = [];
end

CUTOFF_DEVS = 2; % removes edges this many standard deviations below average
TOL = 0.99; % prevent errors by ensuring Hmin <= TOL*Hmax

%% parse/load file inputs
%===========================================================
if isstring(geom) || ischar(geom)
	% check input corresponds to a valid 
	geom = string(geom);
	if ~isscalar(geom)
		error("ERROR: text inputs must correspond to a singular file-name.")
	elseif (exist(geom,"file") ~= 2)
		error("ERROR: the input does not correspond to a file on the current path.")
	end
	
	[~,~,extension] = fileparts(geom);
	extension = lower(extension);
	switch extension
		case ".stl"
			% load the .STL file as a surface triangulation
			geom = stlread(geom);
		case ".step"
			% convert the .STEP file to a fegeometry object
			geom = fegeometry(geom);
		otherwise
			error("ERROR: the input must be an STL or STEP file.")
	end
else
	extension = "none";
end


%% set default parameters
%===========================================================

if isempty(opt.normalize)
	% set domain normalization/centering default
	switch extension
		case ".stl";  opt.normalize = true;
		case ".step"; opt.normalize = false;
		case "none";  opt.normalize = false;
	end
end

if isempty(opt.cullEdges)
	% toggle badly-scaled edge removal default
	switch extension
		case ".stl";  opt.cullEdges = true;
		case ".step"; opt.cullEdges = false;
		case "none";  opt.cullEdges = false;
	end
end

if isempty(opt.maxIters)
  % use fewer refinement iterations by default for STL/STEP files
	switch extension
		case ".stl";  opt.maxIters = 10;
		case ".step"; opt.maxIters = 10;
		case "none";  opt.maxIters = 20;
	end
end

%% remove bad edges from triangulations
%===========================================================
if opt.cullEdges && matches(class(geom),["triangulation","delaunayTriangulation"])
	% get the length of every edge in the mesh
	pts = geom.Points;
	faces = geom.ConnectivityList;
	allEdges = subfaces(faces.',2);
	edgeLen = hypervolumes(pts.',allEdges.');
	
	% determine cutoff for the allowed edge-size
	Hmax = max(edgeLen);
	edgeDist = sqrt(edgeLen/Hmax);
	cutoff = mean(edgeDist) - CUTOFF_DEVS * std(edgeDist);
	cutoff = max(cutoff,0);
	cutoff = Hmax * cutoff^2;
	
	% remove any extremely small edges
	[Hmin,indx] = min(edgeLen);
	while (Hmin < cutoff)
		% replace smallest edge with a single node at the midpoint
		edge_ii = allEdges(indx,:);
		onEdge = ismember(faces,edge_ii);
		
		pts(end+1,:) = mean(pts(edge_ii,:),1); %#ok<AGROW>
		faces(onEdge) = size(pts,1);
		faces(sum(onEdge,2)>1,:) = [];
		
		keepIndex = setdiff(1:size(pts,1),edge_ii);
		[pts,faces] = mesh_subset(pts.',faces.',"pts",keepIndex);
		pts = pts.'; faces = faces.';
		
		% update Hmin
		allEdges = subfaces(faces.',2);
		edgeLen = hypervolumes(pts.',allEdges.');
		[Hmin,indx] = min(edgeLen);
	end

	% re-construct the triangulation
	geom = triangulation(faces,pts);
end


%% convert input to a fegeometry object and determine Hmin
%===========================================================
if ~isa(geom,"fegeometry")
	% convert geometry to a fegeometry object
	geom = fegeometry(geom);
end

% determine Hmin
geomRep = geom.GeomRep;
if isa(geomRep,"pde.DiscreteGeometry")
	% find minimum edge length
	pts = geomRep.Triangulation.Points.';
	tri = geomRep.Triangulation.Triangles.';
	Hmin = min(hypervolumes(pts,subfaces(tri,2).'));
elseif isa(geomRep,"pde.AnalyticGeometry")
	% find minimum distance between any existing vertices
	pts = geomRep.Vertices.';
	[ii,jj] = ncombsk(size(pts,2),2);
	Hmin = min(nonzeros(vecnorm(pts(:,ii) - pts(:,jj),2,1)));
else
	% cannot determine Hmin automatically
	Hmin = 0;
end

% ensure shortest edge in geometry is included
Hmin = max(0,Hmin-2*eps(Hmin));

%% generate the minimally-refined mesh
%===========================================================
Hmax = realmax;
mesh = generateMesh(geom,"Hmax",Hmax,"Hmin",max(0,min(Hmin,TOL*Hmax)), ...
	"Hgrad",opt.Hgrad,"GeometricOrder","linear");
mesh = mesh.Mesh;
pts = mesh.Nodes;
tri = mesh.Elements;
currSize = size(pts,2);


%% update Hmax to improve agreement with the requested mesh size
%===========================================================
% exit if the minimally-refined mesh exceeds the requested size
if (currSize >= nodes) || (opt.maxIters==0)
	obj = FE_mesh(pts,tri);
	return
else
	% determine the current node spacing
	HmaxOld = max(hypervolumes(pts,subfaces(tri,2).'));
end

% determine the ideal node spacing
numDim = size(pts,1);
if numDim~=2 && numDim~=3
	error("ERROR: geometry has an invalid number of dimensions.")
elseif numDim == 2
	area = sum(hypervolumes(pts,tri));
	perim = sum(hypervolumes(pts,boundary_faces(tri).'));
	HmaxNew = node_spacing(nodes,[],area,perim);
elseif numDim == 3
	vol = sum(hypervolumes(pts,tri));
	area = sum(hypervolumes(pts,boundary_faces(tri).'));
	HmaxNew = node_spacing(nodes,vol,area,[]);
end
Hmax = (HmaxNew + HmaxOld)/2;

% update node-separation until the mesh size is optimized
%-----------------------------------------------------------
beta = 1.50; % scaling factor for updating Hmax
thresh = ceil(opt.tolerance * nodes); % convergence threshold
bestSize = currSize;
bestMesh = mesh;

for iter = 1:opt.maxIters
	% generate a mesh
	mesh = generateMesh(geom,"Hmax",Hmax,"Hmin",max(0,min(Hmin,TOL*Hmax)), ...
		"Hgrad",opt.Hgrad,"GeometricOrder","linear");
	mesh = mesh.Mesh;

	prevSize = currSize;
	currSize = size(mesh.Nodes,2);

	if abs(nodes-currSize) < abs(nodes-bestSize)
		% update best result
		bestSize = currSize;
		bestMesh = mesh;
	end

	% update Hmax
	if currSize < (nodes - thresh)
		% Hmax too big
		if prevSize > nodes; beta = sqrt(beta); end
		Hmax = Hmax/beta;
	elseif currSize > (nodes + thresh)
		% Hmax too small
		if prevSize < nodes; beta = sqrt(beta); end
		Hmax = Hmax * beta;
	else
		% exit if mesh-size is within tolerance
		break
	end

	% exit if updates are too small
	if (beta <= 1+opt.converged) 
		break
	end
end

%% return FE_mesh
%===========================================================
% use the mesh where the number of points is closest to the requested value
pts = bestMesh.Nodes;
tri = bestMesh.Elements;

if opt.normalize
	% normalize and center the mesh
	ptsMin = min(pts,[],2);
	ptsMax = max(pts,[],2);
	ptsMid = (ptsMin+ptsMax)/2;
	ptsScl = 1/max(ptsMax-ptsMin);
	pts = ptsScl * (pts - ptsMid);
end
obj = FE_mesh(pts,tri);

%% DONE
end


