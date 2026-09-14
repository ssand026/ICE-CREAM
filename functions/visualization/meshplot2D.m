function [varargout] = meshplot2D(pts,tri,data,opt)
% MESHPLOT2D Visualize multiple scalar fields on the same 2D triangular mesh.
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
% The optional outputs FIG and PLT correspond to the generated figure handle, the
% boundary-mesh patch object handle, and the data plot's graphics object respectively. 
%
% Syntax
% =======
% MESHPLOT2D(PTS,TRI,DATA)
% MESHPLOT2D(___,Name,Value)
% [FIG,PLT] = MESHPLOT2D(PTS,TRI,DATA)
% [FIG,PLT] = MESHPLOT2D(___,Name,Value)
%
% Input Arguments
% ================
% PTS  - (2 x numPts) array containing the coordinates of each mesh vertex
% TRI  - (3 x numTri) array containing the node indices forming each tetrahedron
% DATA - Numeric array corresponding to a single/several scalar fields. Valid formats include:
%		   (3 x numTri x numData) : plots multiple discontinuous scalar fields
%          (numPts x numData)     : plots multiple continuous scalar fields
%          empty array            : plot only the mesh boundary
%
% Name-Value Arguments
% =====================
% style - The data visualization method.
%            "2D" : display top-down view of data (default)
%            "3D" : display data values along the z-axis
%
% cmapLimits - Two-element array specifying of the colormap limits. Values may be numeric or
%              statistical expressions such as [-1,+1],["MIN","MAX"], or ["AVG-2*STD","AVG+2*STD"].
%              The special value ["?","?"] selects a default based on the values in DATA.
% fixedCLims - Logical flag indicating whether the colormap limits are evaluated independently for
%              each data frame.  
% meshStyle - Cell array that passes additional name-value options to the boundary PATCH object.
% axesStyle - Controls the appearance of the coordinate axes:
%	             "default" : standard axes with x, and y labels.
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
% PLT - handle to the patch object containing the displayed data
%
%
% SEE ALSO: MESHPLOT3D, PATCH, STATS_FUNC
arguments
	pts (2,:) double {mustBeReal,mustBeFinite}
	tri (:,:) double {mustBeInteger,mustBePositive}
	data (:,:,:,:) double
	opt.style (1,1) {mustBeMember(opt.style,["2D","3D"])} = "2D";
	opt.cmapLimits (1,2) {validateattributes(opt.cmapLimits,{'string','numeric'},{})} = ["MIN","MAX"];
	opt.fixedCLims (1,1) logical = false;
	opt.axesStyle (1,1) {mustBeMember(opt.axesStyle,["default","none"])} = "default";
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
	isVertexData = isequal(size(data,[1,2]), [numVtx,numTri]);
	
	% setup axes for disjointed data
	%-----------------------------------------------------------
	if isVertexData
		% circumvent the costly overhead of "unique(...)" using:
		%	[C,ia,ic] = matlab.internal.math.uniquehelper(A,doSort,isFirst,byRows)
		fast_unique = @(x,opts) matlab.internal.math.uniquehelper(x,opts(1),opts(2),opts(3));

		[vtx,~,ic] = fast_unique([double(tri(:)),reshape(data,[],numFrames)],[false,true,true]);
		pts = pts(:,vtx(:,1));
		tri = reshape(ic,size(tri));
		data = vtx(:,2:end);
	end
	
	%% handle complex inputs
	%===========================================================
	if ~isreal(data)
		opt.style = opt.style + "_complex";
		
		% set the default colormap properties
		opt.colormap = hsv;
		opt.cmapLimits = (pi/2)*[-1,+1];
		
		% process the data into phase/magnitude
		colorData = angle(data); % phase information
		data = abs(data).^2; % magnitude information
	else
		colorData = data;
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

end

%% figure configs
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
ax.InnerPosition = [0.01 0.01 0.98 0.98];
ax.DataAspectRatio = [1 1 1];
ax.PlotBoxAspectRatio = [1 1 1];
ax.Box = "on";
ax.Clipping = "off";
ax.TickDir = "none";

% axis bounding-box setup
pad = 0.05;
[xyMin, xyMax] = bounds(pts,2);
xyMid = (xyMin + xyMax)/2;
xyLen = (xyMax - xyMin)/2;
xyLen = xyLen + eps(xyLen);
xyBox = xyMid + (1+pad)*[-xyLen,+xyLen];

ax.XLim = xyBox(1,:);
ax.YLim = xyBox(2,:);

% configure axes style
%-----------------------------------------------------------
if (opt.axesStyle == "default")
	% add axes labels
	ax.XTick = xyMid(1);  ax.XTickLabel = "x";
	ax.YTick = xyMid(2);  ax.YTickLabel = "y";
	ax.ZTick = [];
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

% evaluate z-axis limits
%-----------------------------------------------------------
if matches(opt.style,["3D","3D_complex"])
	if opt.fixedCLims
		[zMin,zMax] = bounds(data,"all");
	else
		[zMin,zMax] = bounds(data,[1,2]);
	end
	zMid = (zMin(:) + zMax(:))/2;
	zLen = (zMax(:) - zMin(:))/2;
	zLen = zLen + eps(zLen);
	zLims = zeros(numFrames,2) + (zMid+(1+pad)*[-zLen,+zLen]);
	
	% ensure z-axis scaling is reasonable
	ax.DataAspectRatio = [1, 1, 3*zLen/max(xyLen)];

	% set camera position
	ax.View = [-45, 15];
	axis vis3d
end

% set default colormap
%-----------------------------------------------------------
try %#ok
	if isempty(opt.colormap)
		% check whether data diverges
		isDivergent = any(sign(CLims(:,1))==-1) && any(sign(CLims(:,2))==+1);
		if isDivergent
			opt.colormap = cmap_dark; % use divergent colormap
		else
			opt.colormap = flipud(cmap_mono); % use monotonic colormap
		end
	end
	fig.colormap = opt.colormap;
	ax.colormap = opt.colormap;
end

%% plot the mesh
%===========================================================
if isempty(data)
	% plot boundary mesh without data
	%-----------------------------------------------------------
	defaultOpts = {"FaceColor","#9a8566","FaceAlpha",1,"EdgeColor","#000000","EdgeAlpha",1};
	plt = patch("Vertices",pts.',"Faces",tri.',"Parent",ax, defaultOpts{:}, opt.meshStyle{:});
else
	% plot data and boundary mesh
	%-----------------------------------------------------------
	% initialize default mesh options
	defaultOpts = {"FaceColor","interp","EdgeColor","#242424","EdgeAlpha",1};

	% initialize plot
	plt = patch("Vertices",pts.',"Faces",tri.',"Parent",ax, defaultOpts{:}, opt.meshStyle{:});
	
	% initialize index and define update functions
	INDEX = 0; setIndex(1);
	fig.UserData.setIndex = @setIndex;
	fig.UserData.getIndex = @getIndex;
end

% display figure and set outputs
fig.Visible = "on";
if (nargout > 0)
	varargout{1} = fig;
	varargout{2} = plt;
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
		case "2D"
			set(plt,"Cdata",colorData(:,INDEX));
		case "2D_complex"
			set(plt,"FaceVertexAlphaData",data(:,INDEX));
			set(plt,"Cdata",colorData(:,INDEX));
		case {"3D","3D_complex"}
			set(ax,"ZLim",zLims(INDEX,:));
			set(plt,"Vertices",[pts.',data(:,INDEX)]);
			set(plt,"Cdata",colorData(:,INDEX));
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