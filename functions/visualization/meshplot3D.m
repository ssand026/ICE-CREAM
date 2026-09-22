function [varargout] = meshplot3D(pts,tri,data,opt)
% MESHPLOT3D Visualize multiple scalar fields on the same 3D tetrahedral mesh.
%
% Creates an interactive figure corresponding to a mesh defined via the node coordinates PTS and
% the mesh connectivity TRI. The numeric array DATA corresponds to the nodal-values of one or more 
% scalar fields.
%
% The currently displayed scalar field (or data frame) can be changed via the mouse scroll 
% wheel (while hovering over the figure) or by via the setIndex function located in the figure 
% handle's UserData, ex. "FIG.UserData.setIndex(3)" will display the third data field from the list
% of node values "data". Similarly the index corresponding to the currently displayed data axis can
% be queried via the command "INDEX = FIG.UserData.getIndex()".
%
% The optional outputs FIG, BND, and PLT correspond to the generated figure handle, the
% boundary-mesh patch object handle, and the data plot's graphics object respectively. 
%
% Syntax
% =======
% MESHPLOT3D(PTS,TRI,DATA)
% MESHPLOT3D(___,Name,Value)
% [FIG,BND,PLT] = MESHPLOT3D(PTS,TRI,DATA)
% [FIG,BND,PLT] = MESHPLOT3D(___,Name,Value)
%
% Input Arguments
% ================
% PTS  - (3 x numPts) array containing the coordinates of each mesh vertex
% TRI  - (4 x numTri) array containing the node indices forming each tetrahedron
% DATA - Numeric array corresponding to a single/several scalar fields. Valid formats include:
%		   (4 x numTri x numData) : plots multiple discontinuous scalar fields
%          (numPts x numData)     : plots multiple continuous scalar fields
%          empty array            : plot only the mesh boundary
%
% Name-Value Arguments
% =====================
% style - The data visualization method.
%            "surface"  : interpolated isosurfaces based on equal data-value intervals (default)
%             "volume"  : interpolated isosurfaces based on equal data-volumes
%             "scatter" : scatter plot displaying the data values at each mesh vertex
%             "segment" : displays data values along each each mesh edge
%             "complex" : similar to segment, but uses color to convey the field's phase and opacity
%                         to convey the field's amplitude, used automatically if DATA is complex.
%
% cmapLimits - Two-element array specifying of the colormap limits. Values may be numeric or
%              statistical expressions such as [-1,+1],["MIN","MAX"], or ["AVG-2*STD","AVG+2*STD"].
%              The special value ["?","?"] selects a default based on the values in DATA.
% surfValues - Display these explicit isosurface values. May be numeric or statistical expressions
%              in a similar manner to "cmapLimits". Setting it to an empty array will determine the
%              isosurface values automatically based if "style" equals "surface" or "volume".
% numSurface - Positive integer specifying the number of automatically drawn isosurfaces. Default is 6.
% faceReduce - Attempts reduction of the isosurface face count by this fraction. Default is 0.
% fixedCLims - Logical flag indicating whether the colormap limits and isosurface values are
%              evaluated independently for each data frame. 
% meshStyle - Cell array that passes additional name-value options to the boundary PATCH object.
% axesStyle - Controls the appearance of the coordinate axes:
%	             "default" : standard axes with x, y, and z labels.
%	             "arrows"  : coordinate arrows and labels.
%                "none"    : axes, ticks, and coordinate labels are hidden.
%
% colormap - Optional colormap specification. If empty, the colormap is selected automatically
%            from a divergent/linear colormap based on the data range. 
% wrapScroll - Logical flag controlling whether the scroll wheel will wrap around to the final
%              data frames or will stop at the first/last entry.
%
% Output Arguments
% =================
% FIG - Figure handle for the generated plot
% BND - handle to the patch object displaying the boundary mesh
% PLT - handle to the graphics object displaying the data plot
%
%
% SEE ALSO: MESHPLOT2D, PATCH, MESH_ISOSURFACE, STATS_FUNC, DATA_CONTOURS 
arguments
	pts (3,:) double {mustBeReal,mustBeFinite}
	tri (:,:) double {mustBeInteger,mustBePositive}
	data (:,:,:,:) double
	opt.style (1,1) {mustBeMember(opt.style,["surface","volume","scatter","segment","complex"])} = "surface";
	opt.cmapLimits (1,2) {validateattributes(opt.cmapLimits,{'string','numeric'},{})} = ["MIN","MAX"];
	opt.surfValues (1,:) {validateattributes(opt.surfValues,{'string','numeric'},{})} = [];
    opt.fixedCLims (1,1) logical = false;
    opt.numSurface (1,1) double {mustBePositive,mustBeInteger} = 6;
	opt.faceReduce (1,1) {mustBeInRange(opt.faceReduce,0,1)} = 0;
    opt.axesStyle (1,1) {mustBeMember(opt.axesStyle,["default","arrows","none"])} = "default";
	opt.meshStyle (1,:) cell = {};
	opt.colormap (:,3) double = [];
    opt.wrapScroll (1,1) logical = true;
	
end

% parse the plotting data
%===========================================================
if isempty(data)
    numFrames = 1;
else
	%% check/prepare nodeData
	%===========================================================
	% ensure data has the correct number of dimensions
	if (ndims(data) > 3)
		% try squeezing data
		data = squeeze(data);
		if (ndims(data) > 3)
			error("ERROR: the given data has too many dimensions.");
		end
	end
	
	% re-arrange data so that the data-frames lie along the third axis
	numPts = size(pts,2);
	numVtx = size(tri,1);
	numTri = size(tri,2);
	
	sz = size(data,[1,2,3]);
	isPts = find(sz == numPts);
	isVtx = find(sz == numVtx);
	isTri = find(sz == numTri);
	isOne = find(sz == 1);
	
	if isequal(sz([1,2]),[numPts,1]) || isequal(sz([1,2]),[numVtx,numTri])
		% data is already in a valid format
		dimOrder = [1,2,3];
	elseif isequal(sz([1,2]),[1,numPts]) || isequal(sz([1,2]),[numTri,numVtx])
		% data requires a permutation of the first two dimensions
		dimOrder = [2,1,3];
	elseif isscalar(isVtx) && isscalar(isTri)
		% rearrange data into pagearray of vertex data
		dimOrder = [isVtx,isTri,setdiff(1:3,[isVtx,isTri])];
	elseif isscalar(isPts) && ~isempty(isOne)
		% rearrange data into pagearray of node data
		dimOrder = [isPts,isOne,setdiff(1:3,[isPts,isOne])];
	else
		error("ERROR: the given data has an invalid size.");
	end
	
	% apply correct ordering
	data = permute(data,dimOrder);
	numFrames = size(data,3);
	
	%% handle complex inputs
	%===========================================================
	if ~isreal(data)
		if (opt.style ~= "complex")
			warning("Data has non-real compoents. Switching to the 'complex' display-style.")
		end
		opt.style = "complex";
		
		% set the default colormap properties
		opt.colormap = hsv;
		opt.cmapLimits = (pi/2)*[-1,+1];
		
	elseif (opt.style == "complex")
		warning("Data is non-complex. Switching to the 'scatter' display-style.")
		opt.style = "scatter";
	end
	
	%% evaluate colormap limits
	%===========================================================
	% set the default colormap limit method
	if isequal(opt.cmapLimits,["?","?"])
		if any(sign(data(:))==+1) && any(sign(data(:))==-1)
			opt.cmapLimits = ["AVG-3*STD","AVG+3*STD"];
		else
			opt.cmapLimits = ["MIN","MAX"];
		end
	end
	
	if isstring(opt.cmapLimits)
		% evaluate data statistics ...
		if opt.fixedCLims
			% globally
			statsEvaluator = stats_func(data,1:3);
		else
			% per data index
			statsEvaluator = stats_func(data,1:2);
		end
		
		% evaluate CLims
		CLims = zeros(numFrames,2);
		CLims(:,1) = statsEvaluator(opt.cmapLimits(1));
		CLims(:,2) = statsEvaluator(opt.cmapLimits(2));
	else
		% cmapLimits is numeric
		CLims = zeros(numFrames,2) + opt.cmapLimits;
	end
	% ensure colormap limits are ascending (cLims(:,1) < cLims(:,2))
	issame = (CLims(:,1)==CLims(:,2));
	isdecr = ~(CLims(:,1)<=CLims(:,2));
	CLims(issame,:) = CLims(issame,:) + 10*eps(CLims(issame,:)).*[-1,+1];
	CLims(isdecr,:) = [CLims(isdecr,2),CLims(isdecr,1)];
	
	
	%% evaluate surface values
	%===========================================================
	if isempty(opt.surfValues)
		% surfValues is empty
		surfValues = [];
	elseif isstring(opt.surfValues)
		% surface values are expressions
		if ~exist("statsEvalFunc","var")
			% evaluate data statistics ...
			if opt.fixedCLims
				% globally
				statsEvaluator = stats_func(data,1:3);
			else
				% per data index
				statsEvaluator = stats_func(data,1:2);
			end
		end
		
		% evaluate surfValues
		if opt.fixedCLims
			surfValues = arrayfun(statsEvaluator,opt.surfValues);
		else
			surfValues = zeros(numFrames,numel(opt.surfValues));
			for ii = 1:size(surfValues,2)
				surfValues(:,ii) = statsEvaluator(opt.surfValues(ii));
			end
		end
	else
		% surface values are numeric
		surfValues = opt.surfValues;
		if opt.fixedCLims
			surfValues = repmat(surfValues,numFrames,1);
		end
	end
	surfValues = sort(surfValues,2);
end

%% set figure configs
%===========================================================
fig = figure("Visible","off");
fig.WindowStyle = "normal";
fig.MenuBar = "none";
if numFrames > 1
	fig.WindowScrollWheelFcn = @scroll_update;
end

% axis settings
%-----------------------------------------------------------
ax = axes;
ax.NextPlot = "replaceChildren";
ax.Units = "normalized";
ax.PositionConstraint = "OuterPosition";
ax.OuterPosition = [0.01 0.01 0.98 0.98];
ax.DataAspectRatio = [1 1 1];
ax.PlotBoxAspectRatio = [1 1 1];
ax.Box = "on";
ax.Clipping = "off";
ax.TickDir = "none";

% axis bounding-box setup
pad = 0.05;
[xyzMin, xyzMax] = bounds(pts,2);
xyzMid = (xyzMin + xyzMax)/2;
xyzLen = (xyzMax - xyzMin)/2;
xyzLen = xyzLen + eps(xyzLen);
xyzBox = xyzMid + (1+pad)*[-xyzLen,+xyzLen];

ax.XLim = xyzBox(1,:);
ax.YLim = xyzBox(2,:);
ax.ZLim = xyzBox(3,:);

% configure axes style
%-----------------------------------------------------------
if (opt.axesStyle == "default")
	% add axes labels
	ax.XTick = xyzMid(1);  ax.XTickLabel = "x";
	ax.YTick = xyzMid(2);  ax.YTickLabel = "y";
	ax.ZTick = xyzMid(3);  ax.ZTickLabel = "z";
else
	% remove axes labels
	ax.XTick = [];
	ax.YTick = [];
	ax.ZTick = [];

	% hide axes box
	ax.Color = "none";
	ax.XColor = "none";
	ax.YColor = "none";
	ax.ZColor = "none";
end

% add axes arrows
if (opt.axesStyle == "arrows")
	hold on;
	uvw = diff(xyzBox,1,2);
	xyz = {xyzBox(1,1), xyzBox(2,1), xyzBox(3,1)};
	quiverOpt = {"LineWidth",2*(ax.LineWidth),"ShowArrowHead","off"};
	quiver3(xyz{:}, uvw(1),0,0, "off", "Color", "r", quiverOpt{:})
	quiver3(xyz{:}, 0,uvw(2),0, "off", "Color", "g", quiverOpt{:})
	quiver3(xyz{:}, 0,0,uvw(3), "off", "Color", "b", quiverOpt{:})
	
	textOpt = {"FontWeight","bold","HorizontalAlignment","center"};
	text((1+pad)*xyzBox(1,2), xyzBox(2,1), xyzBox(3,1), "x", "Color", 0.7*[1 0 0], textOpt{:},"VerticalAlignment","baseline")
	text(xyzBox(1,1), (1+pad)*xyzBox(2,2), xyzBox(3,1), "y", "Color", 0.7*[0 1 0], textOpt{:},"VerticalAlignment","middle")
	text(xyzBox(1,1), xyzBox(2,1), (1+pad)*xyzBox(3,2), "z", "Color", 0.7*[0 0 1], textOpt{:},"VerticalAlignment","baseline")
	hold off;
end

% set camera position
%-----------------------------------------------------------
ax.View = [-45, +15];
axis vis3d

% set default colormap
%-----------------------------------------------------------
try %#ok
	if isempty(opt.colormap)
		% check whether data diverges
		isDivergent = any(sign(CLims(:,1))==-1) && any(sign(CLims(:,2))==+1);
		if isDivergent
			opt.colormap = cmap_light; % use divergent colormap
		else
			opt.colormap = flipud(cmap_mono); % use monotonic colormap
		end
	end
	fig.Colormap = opt.colormap;
	ax.Colormap = opt.colormap;
end

%% plot the boundary mesh
%===========================================================
[external] = boundary_faces(tri);
if isempty(data)
	% plot the external boundary and then exit
	
	% generate wireframe boundary mesh
	defaultOpts = {"FaceColor","#9a8566","EdgeColor","#000000","FaceAlpha",0.75,"EdgeAlpha",1};
	bnd = patch("Vertices",pts',"Faces",external,"Parent",ax, defaultOpts{:}, opt.meshStyle{:});
	
	% return the figure handle and exit
	fig.Visible = "on";
	if (nargout > 0)
		varargout{1} = fig;
        varargout{2} = bnd;
        varargout{3} = [];
	end
	return
else
	% generate a transparent boundary mesh
	defaultOpts = {"FaceColor","#242424","EdgeColor","none","FaceAlpha",0.03,"FaceOffsetBias",+1e4};
	bnd = patch("Vertices",pts',"Faces",external,"Parent",ax, defaultOpts{:}, opt.meshStyle{:});
end


%% interpolating contour plots
%===========================================================
if matches(opt.style,["surface","volume"])
	surfPts = cell(1,numFrames);
	surfFace = cell(1,numFrames);
	colorData = cell(1,numFrames);
	
	for jj = 1:numFrames
		% get data/surfaces for the current frame
		data_jj = data(:,:,jj);

		if isempty(surfValues)
			surf_jj = [];
		elseif size(surfValues,1) == 1
			surf_jj = surfValues;
		elseif size(surfValues,1) == numFrames
			surf_jj = surfValues(jj,:);
		end
		
		% determine the ideal contour-surface values for each frame
		%***********************************************************
		if ~isempty(surfValues)
			% use the explicitly given surface values
			svals = surf_jj;
		elseif (opt.style == "surface")
			% contours enclose equal data ranges
			svals = data_contours(data_jj,opt.numSurface,"values");
		elseif (opt.style == "volume")
			% contours enclose equal fractions of the input data
			svals = data_contours(data_jj,opt.numSurface,"volume");
		end
		
		% find the isosurface for each contour value
		%***********************************************************
		if ~isempty(svals)
			svals = sort(svals(:),1,"ascend","ComparisonMethod","abs");
			[newPts,newFace] = mesh_isosurfaces(pts,tri,data_jj,svals);
			
			if (opt.faceReduce ~= 0)
				% reduce the number of faces in each contour
				for kk = 1:numel(svals)
					[newFace{kk},newPts{kk}] = reducepatch(newFace{kk},newPts{kk},1-opt.faceReduce);
				end
			end
			
			% combine contours
			for kk = 1:numel(svals)
				surfFace{jj} = [surfFace{jj}; newFace{kk}+size(surfPts{jj},1)];
				surfPts{jj} = [surfPts{jj}; newPts{kk}];
				colorData{jj} = [colorData{jj}; kk * ones(size(newFace{kk},1),1)];
			end
			colorData{jj} = svals(colorData{jj});
		end
	end
	
	% initialize patch
	%***********************************************************
	plt = patch("Parent",ax,"FaceColor","flat","EdgeColor","none","FaceAlpha",0.3);
	% ax.SortMethod = "childorder";
end


%% scatter/points style plot
%===========================================================
if opt.style == "scatter"
	
	% allow display of scatter plot
	ax.NextPlot = "add";
	
	% plot points/dots
	plt = scatter3(pts(1,:).',pts(2,:).',pts(3,:).',"Parent",ax, ...
		"Marker","o","MarkerFaceColor","flat","MarkerFaceAlpha",1, ...
		"LineWidth",1,"MarkerEdgeAlpha",0);
	
	ax.NextPlot = "replacechildren";
	
	% normalize data to get dot sizes
	dotScale = 32; % max dot size (in pt)
	dataMed = median(data,[1,2]);
	dataStd = std(data,0,[1,2]) + eps;
	sizeData = abs(data - dataMed)./dataStd;
	sizeData = sqrt(min(sizeData,4));
	sizeData = sizeData./max(sizeData,[],[1,2]);
	sizeData = floor(dotScale*sizeData);
	sizeData(sizeData <= 0) = eps;
	
	% set color map
	colorData = data;
end

%% segment/wireframe style plot
%===========================================================
if opt.style == "segment"
	[allFaces] = subfaces(tri,3);
	[internal] = setdiff(allFaces,external,"rows");
	
	% plot mesh edges
	plt = patch("Vertices",pts.',"Faces",internal,"Parent",ax, ...
		"EdgeColor","interp","EdgeAlpha","flat","lineWidth",3, ...
		"FaceColor","none","FaceAlpha",0);
	
	colorData = data; % set color map
	alphaData = abs(data); % set alpha map
	
end

%% complex data plot
%===========================================================
if opt.style == "complex"
	[allFaces] = subfaces(tri,3);
	[internal] = setdiff(allFaces,external,"rows");
	
	% plot mesh edges
	plt = patch("Vertices",pts.',"Faces",internal,"Parent",ax, ...
		"EdgeColor","interp","EdgeAlpha","flat","lineWidth",3, ...
		"FaceColor","none","FaceAlpha",0);
	
	colorData = angle(data); % phase information
	alphaData = abs(data).^2; % magnitude information
end


%% display figure and define update functions
%===========================================================
INDEX = 0; setIndex(1);
fig.UserData.setIndex = @setIndex;
fig.UserData.getIndex = @getIndex;
fig.Visible = "on";
if (nargout > 0)
	varargout{1} = fig;
	varargout{2} = bnd;
    varargout{3} = plt;
end

%% internal functions
%***********************************************************
	function scroll_update(~,event)
	% use mouse scroll wheel to update the displayed data index
	newIndex = INDEX + (event.VerticalScrollCount);
	setIndex(newIndex)
	end
%***********************************************************
	function [out] = getIndex()
	% return the current INDEX value and reset UserData.index
	fig.UserData.index = INDEX;
	out = INDEX;
	end
%***********************************************************
	function setIndex(newIndex)
	%% check if the 'newIndex' value is valid
	%-----------------------------------------------------------
	try; newIndex = double(newIndex); end %#ok
	if isnumeric(newIndex) && isscalar(newIndex) && isfinite(newIndex)
		% convert newIndex to an integer value between [1,numFrames]
		newIndex = round(double(newIndex));
		if opt.wrapScroll
			newIndex = mod(newIndex,numFrames) + (newIndex < 0);
			newIndex(newIndex==0) = numFrames;
		else
			newIndex = min(max(1,newIndex),numFrames);
		end
	else
		% newIndex is invalid, exit without updating INDEX
		return
	end
	
	% check if INDEX has changed
	if (INDEX == newIndex)
		% exit without updating
		return
	else
		INDEX = newIndex;
		fig.UserData.index = INDEX;
	end
	
	%% update the figure
	%-----------------------------------------------------------
	set(plt,"AlphaDataMapping","scaled");
	switch opt.style
		case {"surface","volume"}
			set(plt,"Vertices",surfPts{INDEX});
			set(plt,"Faces",surfFace{INDEX});
			set(plt,"CData",colorData{INDEX});
		case "scatter"
			set(plt,"SizeData",sizeData(:,INDEX));
			set(plt,"CData",colorData(:,INDEX));
		case {"segment","complex"}
			set(plt,"FaceVertexAlphaData",alphaData(:,INDEX));
			set(plt,"CData",colorData(:,INDEX))
	end
	
	% update color-data limits
	set(ax,"CLim",CLims(INDEX,:));
	
	% update figure title
	if numFrames > 1
		figName = sprintf("  [%i / %i]",INDEX,numFrames);
		set(fig,"Name",figName);
	end
	%-----------------------------------------------------------
	end
%***********************************************************
%% DONE
end