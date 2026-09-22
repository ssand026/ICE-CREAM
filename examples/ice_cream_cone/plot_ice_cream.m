%% load system
clc; clearvars;
load("ice_cream_QC.mat")


%% plot the ice cream mesh
%-----------------------------------------------------------
fig = meshplot3D(mesh.pts,mesh.tri,[]);
axOrig = fig.CurrentAxes;
axOrig.CameraViewAngleMode = "auto";
bbox = boundingBox(axOrig);
boxW = bbox(3); boxH = bbox(4);
fig.Position = fig.Position .* [1, 1, sqrt(boxW/boxH), sqrt(boxH/boxW)];

% get normalized z-positions/ mesh regions
z = mesh.pts(3,:);
z = (z(:)-min(z))/(max(z)-min(z));
z1 = 0.401; % 1st scoop height
z2 = 0.668; % 2nd scoop height
reg = 1 + (z > z1) + (z > z2);

% color mesh to look like ice-cream
plt = fig.CurrentAxes.Children;
plt.FaceColor = 'interp';
plt.FaceAlpha = 1;

plt.FaceVertexCData = reg(:);
colormap([ ...
	0.53, 0.40, 0.25; ... % brown
	0.90, 0.40, 0.40; ... % pink
	0.95, 0.95 ,0.80 ... % white
	])

% SAVE
drawnow;
exportgraphics(fig,"ice_cream_mesh.pdf","Padding","figure","Resolution",500);

%% eigenstates plot
%-----------------------------------------------------------
clc;
close all;

% generate eigenstate plots
mpts = mesh.pts;
mtri = mesh.tri;
fig = meshplot3D(mpts,mtri,to_basis(states,"full",BCM{:}),"style","surface", ...
	"cmapLimits",["-3*STD","+3*STD"],"fixedClims",0,"axesStyle","none","numSurface",8);
axOrig = fig.CurrentAxes;
axOrig.CameraViewAngleMode = "auto";
bbox = boundingBox(axOrig);

% generate eigenstate grid
nrow = 2;
ncol = 7;

boxW = bbox(3); boxH = bbox(4);
xyAspect = (ncol * boxW)/(nrow * boxH);

grid = figure;
grid.MenuBar = "none";
grid.Position = grid.Position .* [1, 1, sqrt(xyAspect), sqrt(1/xyAspect)];

for ii = 1:nrow
	for jj = 1:ncol
		% copy relevant axes to grid figure
		indx = jj + (ii-1) * ncol;
		fig.UserData.setIndex(indx)
		ax = copyobj(fig.CurrentAxes,grid);
		
		% set axes settings
		ax.Box = "off";
		ax.XAxis.Visible = "off";
		ax.YAxis.Visible = "off";
		ax.ZAxis.Visible = "off";
		ax.Units = "normalized";
		
		% reposition and label
		ax.OuterPosition = [(jj-1)/ncol, 0.02+(nrow-ii)/nrow, 1/ncol, 1/nrow];
		
		ax.Title.String = "E_{"+indx+"} = "+sprintf("%1.2f",E(indx)/(1e-3*eV))+" meV";
		ax.Title.Units = "normalized";
		ax.Title.Position = [0.5, -0.025, 0];
	end
end

% SAVE
drawnow;
% exportgraphics(grid,"ice_cream_eigs.pdf","Padding","figure","Resolution",500);


%% compare electric/magnetic TDM
%-----------------------------------------------------------
numEigs = size(states,2);
TDM_E = get_TDM(states(:,1:numEigs),muE,"noDiags",true); TDM_E = vecnorm(TDM_E,2,3);
TDM_B = get_TDM(states(:,1:numEigs),muB,"noDiags",true); TDM_B = vecnorm(TDM_B,2,3);
TDM_E = TDM_E/max(TDM_E(:));
TDM_B = TDM_B/max(TDM_B(:));

% options
ar = [1 0.45]; % aspect ratio w,h
m = 0.00; % margin
nc = 10; % colors per magnitude step
cs = 2; % scaling exponent for interpolation (1 == linear)
cmap = cmap_mono;
cmap = cmap(4:end,:);

% figure
fig = figure;
fig.Units = "points";
fig.Position(3:4) = ar .* fig.Position(3:4);

% electric
%-----------------------------------------------------------
cMag = -4; % colorscale minimum

ax1 = axes;
img = log10(TDM_E); 
imagesc(img); axis image; 
clim([cMag 0])

xline((2:numEigs)-1/2,"LineWidth",1)
yline((2:numEigs)-1/2,"LineWidth",1)
xlabel("$\vert{\Psi_i}\rangle$","Interpreter","latex");
ylabel("$\langle{\Psi_f}\vert$","Interpreter","latex","Rotation",0);
xticks(1:numEigs); xtickangle(0)
yticks(1:numEigs)

ax1.Title.String = "(a)";
ax1.TickDir = "none";
ax1.Units = "normalized";
ax1.OuterPosition = [m, 0, 0.5-2*m, 1];

% colorbar
cbar = colorbar;
cbar.TickLabels = "10^{"+(cMag:0)+"}";
cbar.Ticks = (cMag:0);

ax1.Colormap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,nc*(0-cMag)).^cs);
ax1.PositionConstraint = "InnerPosition";
cbar.Position(3) = 0.02;
ax1.PositionConstraint = "OuterPosition";

% magnetic
%-----------------------------------------------------------
cMag = -4;

ax2 = axes;
img = log10(TDM_B); 
imagesc(img); axis image
clim([cMag 0])

xline((2:numEigs)-1/2,"LineWidth",1)
yline((2:numEigs)-1/2,"LineWidth",1)
xlabel("$\vert{\Psi_i}\rangle$","Interpreter","latex");
ylabel("$\langle{\Psi_f}\vert$","Interpreter","latex","Rotation",0);
xticks(1:numEigs); xtickangle(0)
yticks(1:numEigs)

ax2.Title.String = "(b)";
ax2.TickDir = "none";
ax2.Units = "normalized";
ax2.OuterPosition = [0.5+m, 0, 0.5-2*m, 1];

% colorbar
cbar = colorbar;
cbar.TickLabels = "10^{"+(cMag:0)+"}";
cbar.Ticks = (cMag:0);

ax2.Colormap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,nc*(0-cMag)).^cs);
ax2.PositionConstraint = "InnerPosition";
cbar.Position(3) = 0.02;
ax2.PositionConstraint = "OuterPosition";

% SAVE
drawnow;
% exportgraphics(fig,"ice_cream_tdm.pdf","Padding","figure","Resolution",500);

%% compare electric/magnetic XYZ
%-----------------------------------------------------------
% options
ar = [1 0.875]; % aspect ratio w,h
m = 0.002; % margin
cMag = -4; % logscale minimum
nc = 10; % colors per magnitude step
cs = 2; % scaling exponent for interpolation (1 == linear)

% get tdm
numEigs = size(states,2);
TDM_E = get_TDM(states(:,1:numEigs),muE,"noDiags",1); 
TDM_B = get_TDM(states(:,1:numEigs),muB,"noDiags",1);

%-----------------------------------------------------------
% colormap
cmap = cmap_mono;
cmap = cmap(3:end,:);
cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,nc*(0-cMag)).^cs);

% electric tdm scaling
E_mag = max(vecnorm(TDM_E,2,3),[],"all");
E_lims = [cMag,0]+log10(E_mag);
E_tick = [floor(E_lims(1)):ceil(E_lims(2))];

% magnetic tdm scaling
B_mag = max(vecnorm(TDM_B,2,3),[],"all");
B_lims = [cMag,0]+log10(B_mag);
B_tick = [floor(B_lims(1)):ceil(B_lims(2))];
%-----------------------------------------------------------
% axis locations
pos = cell(1,4);
pos{1} = [m, 0.5+m, 0.5-2*m, 0.5-2*m];
pos{2} = [0.5+m, 0.5+m, 0.5-2*m, 0.5-2*m];
pos{3} = [m, m, 0.5-2*m, 0.5-2*m];
pos{4} = [0.5+m, m, 0.5-2*m, 0.5-2*m];

% axis data
imgs = cell(1,4);
imgs{1} = vecnorm(TDM_E(:,:,[1,2]),2,3);
imgs{2} = vecnorm(TDM_E(:,:,[3]),2,3);
imgs{3} = vecnorm(TDM_B(:,:,[1,2]),2,3);
imgs{4} = vecnorm(TDM_B(:,:,[3]),2,3);
imgs = cellfun(@log10,imgs,"UniformOutput",false);

% axis options
titles = {"(a)","(b)","(c)","(d)"};
ticks = {E_tick, E_tick, B_tick, B_tick};
limits = {E_lims, E_lims, B_lims, B_lims};
%-----------------------------------------------------------

% figure
fig = figure;
fig.Units = "points";
fig.Position(3:4) = ar .* fig.Position(3:4);

% generate each axis
for ii = 1:4
	ax = axes;
	imagesc(imgs{ii}); axis image
	ax.YDir = "normal";
	clim(limits{ii})
	
	xline((2:numEigs)-1/2,"LineWidth",1,"Alpha",1)
	yline((2:numEigs)-1/2,"LineWidth",1,"Alpha",1)
	
	xlabel("$\vert{\Psi_\mathrm{i}}\rangle$","Interpreter","latex");
	ylabel("$\langle{\Psi_\mathrm{f}}\vert\quad$","Interpreter","latex", ...
		"Rotation",0,"HorizontalAlignment","center");
	xticks(1:numEigs); xtickangle(0)
	yticks(1:numEigs)
	xlim([1/2,numEigs+1/2])
	ylim([1/2,numEigs+1/2])
	xticklabels("\fontsize{10}"+(1:numEigs))
	yticklabels("\fontsize{10}"+(1:numEigs))
	
	ax.Title.String = titles{ii};
	ax.TickDir = "none";
	ax.Units = "normalized";
	ax.OuterPosition = pos{ii};
	
	% colorbar
	cbar = colorbar;
	cbar.TickLabels = "10^{"+ticks{ii}+"}";
	cbar.Ticks = ticks{ii};
	
	ax.Colormap = cmap;
	ax.PositionConstraint = "InnerPosition";
	cbar.Position(3) = 0.02;
	ax.PositionConstraint = "OuterPosition";
end

% SAVE
pause(0.5)
drawnow;
% exportgraphics(fig,"ice_cream_tdm.pdf","Padding","figure","Resolution",500);


%%
% fig = figure;
% cOrd = colororder("gem12"); 
% cOrd = cOrd([7,10,2,3,8,5,11,6,1,9,4],:);
% % cOrd = cOrd([9,2,6, 8,3,5, 10,1,7, 11, 4],:);
% colororder(fig,cOrd)
% plot((1:size(cOrd,1))/3+sinpi(linspace(0,7,100)).',"LineWidth",12)

%% plot QC magnetic
%-----------------------------------------------------------
% units
[ps,mT] = atomic_units("ps","mT");
fld = field.'/mT;

fig = figure;
cOrd = colororder("gem12"); cOrd = cOrd([1,10,5],:);
colororder(fig,cOrd)
fig.Position = fig.Position .* [1 1 1 0.8];
nRow = 3;
marg = 0.02;

% field
%-----------------------------------------------------------
ax1 = axes;
ax1.Units = "normalized";
ax1.OuterPosition = [0, 2/nRow, 1, (1/nRow-marg/2)];

x = tgrid/ps;
y = fld;
plot(x,y)
xlim([0,x(end)])
ylim(80*[-1,+1])
xticks(0:x(end))
ax1.XGrid = "on";
ax1.YGrid = "on";

xlabel("Time (ps)","VerticalAlignment","Top")
ylabel("Magnetic Field (mT)")
legend({"B{\fontsize{8}x}","B{\fontsize{8}y}","B{\fontsize{8}z}"}, ...
	"IconColumnWidth",24,"FontSize",10,"location","southeast")

% P(t)
%-----------------------------------------------------------
ax2 = axes;
ax2.Units = "normalized";
ax2.OuterPosition = [0, 1/nRow, 1, (1/nRow-marg/2)];

x = tgrid/ps;
y = squeeze(Pt);
plot(x,y)
xlim([0,x(end)])
xticks(0:x(end))
yticks(0:0.2:1)

ax2.XGrid = "on";
ax2.YGrid = "on";
xlabel("Time (ps)")
ylabel("Transition Probability")

% fft
%-----------------------------------------------------------
ax3 = axes;
ax3.Units = "normalized";
ax3.OuterPosition = [0, 0/nRow, 1, (1/nRow-marg/2)];

t = tgrid/atomic_units("s");
[freq,amp] = fft_plot(fld,t,4);
x = freq/1e9;
y = amp/max(amp(:));
plot(x,y.*[1,0.95,1])

xticks(0:200:x(end))
xlim([0,3000])
ylim([0,1]+0.025*[-1,+1])
ax3.XGrid = "on";


xlabel("Frequency (GHz)")
ylabel("Magnitude (arb.)")
legend({"B{\fontsize{8}x}","B{\fontsize{8}y}","B{\fontsize{8}z}"}, ...
	"IconColumnWidth",24,"FontSize",10)

% SAVE
drawnow;
% exportgraphics(fig,"ice_cream_QC.pdf","Padding","figure","Resolution",500);

%% time-dependent density
%-----------------------------------------------------------
close all;

numSlice = round(numTau/20);
slices = round(linspace(1,numel(tgrid),numSlice+1));
slices = slices(2:end);
tSlice = round(tgrid(slices)/fs,1);
rho_t = squeeze(abs(psiFwd(:,:,slices)).^2);


fig = meshplot3D(mesh.pts,mesh.tri,to_basis(rho_t,"full",BCM{:}), ...
	"style","volume","numSurface",8,"axesStyle","none",...
	"cmapLimits",["0","AVG+4*STD"],"fixedClims",true,"faceReduce",0);

view(90,90)
axOrig = fig.CurrentAxes;
axOrig.CameraViewAngleMode = "auto";
bbox = boundingBox(axOrig);

% add geometry outline
% hold on
% indx = all(zEdg <= 0.1,1);
% plot3(xEdg(:,indx),yEdg(:,indx),zEdg(:,indx),"-","LineWidth",1.5,"Color",0.35*[1 1 1])
% hold off


% set patch display options
% axOrig.Children = flipud(axOrig.Children);
axOrig.SortMethod = "childorder";
colormap(cmap_mono)

% set axis background color
pat = axOrig.Children(end);
pat.FaceAlpha = 1;
pat.FaceColor = 0*[1 1 1];

return

% generate timeslice grid
%-----------------------------------------------------------
nrow = 3;
ncol = numSlice/nrow;

boxW = bbox(3); boxH = bbox(4);
xyAspect = (ncol * boxW)/(nrow * boxH);

grid = figure;
grid.MenuBar = "none";
grid.Position = grid.Position .* [1, 1, 1, (1/xyAspect)];

for ii = 1:nrow
	for jj = 1:ncol
		% copy relevant axes to grid figure
		indx = jj + (ii-1) * ncol;
		fig.UserData.setIndex(indx)
		ax = copyobj(fig.CurrentAxes,grid);
		
		% set axes settings
		ax.Box = "off";
		ax.Toolbar.Visible = "off";
		ax.XAxis.Visible = "off";
		ax.YAxis.Visible = "off";
		ax.ZAxis.Visible = "off";
		ax.Units = "normalized";

		% reposition and label
		ax.InnerPosition = [(jj-1)/ncol, (nrow-ii)/nrow, 1/ncol, 1/nrow];
		ax.Title.String = ""+tSlice(indx)+" fs";
		ax.Title.Color = "w";
		ax.Title.Units = "normalized";
		ax.Title.Position = [0.5, 0.02, 0];
		ax.Title.FontWeight = "normal";
	end
end

% SAVE
drawnow;
% exportgraphics(grid,"ice_cream_rho.pdf","Padding","figure","Resolution",500)

%% external functions
%===========================================================
function [freq, amp] = fft_plot(x,t,pad,dim)
% performs the real-space fourier transform of an input array/waveform
% outputs the array of frequencies and amplitudes
arguments
	x (:,:) double
	t (1,:) double
	pad (1,1) double {mustBeInteger,mustBeNonnegative} = 0;
	dim (1,:) double {mustBeInteger,mustBePositive} = [];
end

% intialize data dim
if isempty(dim)
	dim = find(size(x)~=1,1,"first");
end

% T = (max(t) - min(t)); % time duration
% numData = size(x,dim); % number of data points
% sampleFreq = numData/T;
% freq = linspace(0,sampleFreq/2,2^nextpow2(numData));
% amp = nufft(x,t,freq,dim);

numData = size(x,dim); % number of data points
numQuery = 2^(nextpow2(numData)+1); % number of interpolation points
numSamp = (2^pad) * numQuery; % number of sample points
numFreq = numSamp/2 + 1; % number of frequency values

t = (t - min(t)); % shift time values
T = (max(t) - min(t)); % time duration
tq = linspace(0,T,numQuery);
xq = interp1(t,x,tq,"cubic");

amp = fft(xq,numSamp,dim);
amp = abs(amp/numSamp);
amp = part(amp,(1:numFreq),dim);

sampleFreq = numQuery/T;
freq = linspace(0,sampleFreq/2,numFreq);
end
%===========================================================
function [bbox] = boundingBox(ax)
% get the 2D bounding box of a 3D axis

% temporarily set the CameraViewAngleMode
cmode = ax.CameraViewAngleMode;
ax.CameraViewAngleMode = "auto";

% corners of the 3D axes box
[X, Y, Z] = ndgrid(ax.XLim, ax.YLim, ax.ZLim);
corners = [X(:), Y(:), Z(:), ones(8,1)].';

% get camera view and projection parameters
[az, el] = view(ax);
fov = ax.CameraViewAngle;
cam_pos = campos(ax);
cam_target = camtarget(ax);

% compute transformation matrix
viewMat = viewmtx(az, el, fov, cam_pos - cam_target); 
projCoords = viewMat * corners;

% convert homogeneous coordinates to 2D
x2d = projCoords(1,:) ./ projCoords(4,:);
y2d = projCoords(2,:) ./ projCoords(4,:);


% get teh normalized bounding box [x_min, y_min, width, height]
w = max(x2d) - min(x2d);
h = max(y2d) - min(y2d);
d = max(w,h);
bbox = [min(x2d)/d+1/2, min(y2d)/d+1/2, w/d, h/d];

% restore the CameraViewAngleMode
ax.CameraViewAngleMode = cmode;
% DONE
end
%===========================================================