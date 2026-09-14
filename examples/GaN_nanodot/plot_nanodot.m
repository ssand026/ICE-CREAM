%% load system
clc; clearvars;
load("nanodot_4_7_QC.mat")

%% plot the nanodot mesh
%-----------------------------------------------------------

fig = meshplot3D(mesh.pts,mesh.tri,[]);
view(-15,20)

axOrig = fig.CurrentAxes;
axOrig.CameraViewAngleMode = "auto";
bbox = boundingBox(axOrig);
boxW = bbox(3); boxH = bbox(4);
fig.Position(3:4) = fig.Position(3) .* [1, (boxH/boxW)];

plt = fig.CurrentAxes.Children;
plt.FaceAlpha = 1;
plt.Faces = align_normals(plt.Faces);
plt.BackFaceLighting="lit";

hold on
plot3(xEdg,yEdg,zEdg,"-k")
hold off

l = light;
lightangle(-90,15)
material([1 1 0.25 5])

% SAVE
drawnow;
% exportgraphics(fig,"nanodot_mesh.pdf","Padding","figure","Resolution",500)

%% eigenstates plot
%-----------------------------------------------------------
clc;
% close all;

% generate eigenstate plots
mpts = mesh.pts;
mtri = mesh.tri;
fig = meshplot3D(mpts,mtri,to_basis(states,"full",BCM{:}),"style","surface", ...
	"cmapLimits",["-4*STD","+4*STD"],"fixedClims",0,"axesStyle","none","surfCount",8);
view(-7.5,35)
axOrig = fig.CurrentAxes;
axOrig.CameraViewAngleMode = "auto";
bbox = boundingBox(axOrig);
hold on
plot3(xEdg,yEdg,zEdg,"-k","LineWidth",1.25)
hold off

% generate eigenstate grid
%-----------------------------------------------------------
nrow = 3;
ncol = 5;

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
		ax.Toolbar.Visible = "off";
		ax.XAxis.Visible = "off";
		ax.YAxis.Visible = "off";
		ax.ZAxis.Visible = "off";
		ax.Units = "normalized";

		% reposition and label
		ax.InnerPosition = [(jj-1)/ncol, (nrow-ii)/nrow, 1/ncol, 1/nrow];
		% ax.Title.String = "\psi_{"+indx+"}";

		ax.Title.String = "E_{"+indx+"} = "+sprintf("%1.3f",E(indx)/eV)+" eV";
		ax.Title.Units = "normalized";
		ax.Title.Position = [0.5, 0.75, 0];
	end
end

% SAVE
drawnow;
% exportgraphics(grid,"nanodot_eigs.pdf","Padding","figure","Resolution",500)

%% compare electric/magnetic TDM
%-----------------------------------------------------------
numEigs = size(states,2);
numEigs = min(numEigs,15);
TDM_E = vecnorm(get_TDM(states(:,1:numEigs),muE,"noDiags",true),2,3);
TDM_B = vecnorm(get_TDM(states(:,1:numEigs),muB,"noDiags",true),2,3);
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
cMag = -4;
tmp_cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,nc*(0-cMag)).^cs);

ax1 = axes;
img = log10(TDM_E); 
imagesc(img); axis image; 
clim([cMag 0])

xline((2:numEigs)-1/2,"LineWidth",1)
yline((2:numEigs)-1/2,"LineWidth",1)
xlabel("$\vert{\psi_i}\rangle$","Interpreter","latex");
ylabel("$\vert{\psi_f}\rangle$","Interpreter","latex");
xticks(1:numEigs); xtickangle(0)
yticks(1:numEigs)

ax1.Title.String = "(a)";
ax1.TickDir = "none";
ax1.Units = "normalized";
ax1.OuterPosition = [m, 0, 0.5-2*m, 1];
ax1.Colormap = tmp_cmap;

% colorbar
cbar = colorbar;
cbar.TickLabels = "10^{"+(cMag:0)+"}";
cbar.Ticks = (cMag:0);
ax1.PositionConstraint = "InnerPosition";
cbar.Position(3) = 0.02;

% magnetic
%-----------------------------------------------------------
cMag = -4;
tmp_cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,nc*(0-cMag)).^cs);

ax2 = axes;
img = log10(TDM_B); 
imagesc(img); axis image
clim([cMag 0])

xline((2:numEigs)-1/2,"LineWidth",1)
yline((2:numEigs)-1/2,"LineWidth",1)
xlabel("$\vert{\psi_i}\rangle$","Interpreter","latex");
ylabel("$\vert{\psi_f}\rangle$","Interpreter","latex");
xticks(1:numEigs); xtickangle(0)
yticks(1:numEigs)

ax2.Title.String = "(b)";
ax2.TickDir = "none";
ax2.Units = "normalized";
ax2.OuterPosition = [0.5+m, 0, 0.5-2*m, 1];
ax2.Colormap = tmp_cmap;

% colorbar
cbar = colorbar;
cbar.TickLabels = "10^{"+(cMag:0)+"}";
cbar.Ticks = (cMag:0);
ax2.PositionConstraint = "InnerPosition";
cbar.Position(3) = 0.02;

% SAVE
drawnow;
exportgraphics(fig,"nanodot_tdm.pdf","Padding","figure","Resolution",500)

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
numEigs = min(numEigs,15);
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
% exportgraphics(fig,"nanodot_tdm_split.pdf","Padding","figure","Resolution",500)

%% plot QC
%-----------------------------------------------------------
% restore field rotation
u = ["mV","nm"];
E_scl = atomic_units(u(1))/atomic_units(u(2));
E_lbl = (u(1)+"/"+u(2));
% fld = (rotT\field).'/E_scl;
fld = field.' ./E_scl;

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

x = tgrid/fs;
y = fld;
plot(x,y)

xlim([0,x(end)])
ylim([-10,+10])
xticks(0:100:1500)
ax1.XGrid = "on";
ax1.YGrid = "on";

xlabel("Time (fs)")
ylabel("Electric Field ("+E_lbl+")")
legend({"E{\fontsize{8}x}","E{\fontsize{8}y}","E{\fontsize{8}z}"}, ...
	"IconColumnWidth",24,"FontSize",10,"location","southeast")

% P(t)
%-----------------------------------------------------------
ax2 = axes;
ax2.Units = "normalized";
ax2.OuterPosition = [0, 1/nRow, 1, (1/nRow-marg/2)];

x = tgrid/fs;
y = squeeze(Pt);
plot(x,y)

xlim([0,x(end)])
xticks(0:100:x(end))
yticks(0:0.2:1)
ax2.XGrid = "on";
ax2.YGrid = "on";

xlabel("Time (fs)")
ylabel("Transition Probability")


% fft
%-----------------------------------------------------------
ax3 = axes;
ax3.Units = "normalized";
ax3.OuterPosition = [0, 0/nRow, 1, (1/nRow-marg/2)];

t = tgrid/atomic_units("s");
[freq,amp] = fft_plot(fld,t,0);
x = freq/1e12;
y = amp/max(amp(2:end,:),[],"all");
plot(x,y)

xticks(0:2:x(end))
xlim([0,50])
ylim([0,1]+0.02*[-1,+1])
ax3.XGrid = "on";

xlabel("Frequency (THz)")
ylabel("Magnitude (arb.)")
legend({"E{\fontsize{8}x}","E{\fontsize{8}y}","E{\fontsize{8}z}"}, ...
	"IconColumnWidth",24,"FontSize",10)


% SAVE
%-----------------------------------------------------------
drawnow;
% exportgraphics(fig,"nanodot_QC.pdf","Padding","figure","Resolution",500)

%% compare QC and PAT
%-----------------------------------------------------------
% restore field rotation
u = ["mV","nm"];
E_scl = atomic_units(u(1))/atomic_units(u(2));
E_lbl = (u(1)+"/"+u(2));


fig = figure;
cOrd = colororder("gem12"); cOrd = cOrd([1,10,5],:);
colororder(fig,cOrd)
fig.Position = fig.Position .* [1 1 1 0.8];


% load both variables
gamma = 1;
load("nanodot_4_7_PAT.mat","Pt","field")
Pt_PAT = Pt;
fld_PAT = gamma * field.' ./E_scl;

load("nanodot_4_7_QC.mat","Pt","field")
Pt_QC = Pt;
fld_QC = field.' ./E_scl;


nRow = 3;
nCol = 2;
marg = 0.02;

% plot fields
%-----------------------------------------------------------
data = {fld_PAT,fld_QC};
mag = [800,50];
for ii = 1:2
	ax = axes("Units","normalized");
	ax.OuterPosition = pbox(nRow,nCol,1,ii,marg);
	
	x = tgrid/fs;
	plot(x,data{ii},"LineWidth",1.5)
	
	xlim([0,x(end)])
	ylim(mag(ii)*[-1,+1])
	xticks(0:30:300)
	yticks(linspace(-mag(ii),+mag(ii),5))
	ax.XGrid = "on";
	ax.YGrid = "on";
	
	xlabel("Time (fs)")
	ylabel("Electric Field ("+E_lbl+")")
	% legend({"E{\fontsize{8}x}","E{\fontsize{8}y}","E{\fontsize{8}z}"}, ...
	% 	"IconColumnWidth",18,"FontSize",10,"location","southeast","orientation","vertical")
end

% plot fft
%-----------------------------------------------------------
for ii = 1:2
	ax = axes("Units","normalized");
	ax.OuterPosition = pbox(nRow,nCol,2,ii,marg);
	
	t = tgrid/atomic_units("s");
	[freq,amp] = fft_plot(data{ii},t,4);
	x = freq/1e12;
	y = amp/max(amp(2:end,:),[],"all");
	plot(x,y,"LineWidth",1.5)
	
	xlim([0,400])
	ylim([0,1]+0.05*[-1,+1])
	xticks(0:50:x(end))
	yticks(0:1)
	ax.XGrid = "on";
	ax.YGrid = "on";

	xlabel("Frequency (THz)")
	ylabel("Magnitude (arb.)")
	legend({"E{\fontsize{8}x}","E{\fontsize{8}y}","E{\fontsize{8}z}"}, ...
		"IconColumnWidth",18,"FontSize",10)
end

% plot P(t)
%-----------------------------------------------------------
ax = axes("Units","normalized");
ax.OuterPosition = pbox(nRow,nCol,3,1:2,marg);

x = tgrid/fs;
y = [squeeze(Pt_PAT),squeeze(Pt_QC)];
plot(x,y)

xlim([0,x(end)])
xticks(0:30:x(end))
yticks(0:0.2:1)
ax.XGrid = "on";
ax.YGrid = "on";

xlabel("Time (fs)")
ylabel("Transition Probability")
legend({"PAT","QC"},"location","northwest","IconColumnWidth",18,"FontSize",10)

% SAVE
%-----------------------------------------------------------
drawnow;
% exportgraphics(fig,"nanodot_PAT_vs_QC.pdf","Padding","figure","Resolution",500)

%% time-dependent density
%-----------------------------------------------------------
close all;

numSlice = 30;%round(numTau/20);
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
hold on
indx = all(zEdg <= 0.1,1);
plot3(xEdg(:,indx),yEdg(:,indx),zEdg(:,indx),"-","LineWidth",1.5,"Color",0.35*[1 1 1])
hold off


% set patch display options
% axOrig.Children = flipud(axOrig.Children);
axOrig.SortMethod = "childorder";
colormap(cmap_mono)

% set axis background color
pat = axOrig.Children(end);
pat.FaceAlpha = 1;
pat.FaceColor = 0*[1 1 1];


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
		ax.Title.FontSize = 10;
	end
end

% SAVE
drawnow;
% exportgraphics(grid,"nanodot_rho.pdf","Padding","figure","Resolution",500)


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
function [pos] = pbox(numRow,numCol,row,col,marg)
% gives normalized positions for the x and y ranges in a plot grid
arguments
	numRow (1,1) double {mustBeInteger,mustBePositive}
	numCol (1,1) double {mustBeInteger,mustBePositive}
	row (1,:) double {mustBeInteger,mustBePositive}
	col (1,:) double {mustBeInteger,mustBePositive}
	marg (1,:) double {mustBeNonnegative} = 0;
end

if numel(marg) > 2
	error("ERROR: margin must contain 2 elements at most.")
elseif isscalar(marg)
	marg = [marg,marg];
end

if max(row) > numRow
	error("ERROR: row index exceeds number of rows.")
end
if max(col) > numCol
	error("ERROR: column index exceeds number of columns.")
end

% north and south most points
yGrid = linspace(1,0,numRow+1);
n = yGrid(min(row)) - marg(2)/2;
s = yGrid(max(row)+1) + marg(2)/2;

% east and west most points
xGrid = linspace(0,1,numCol+1);
e = xGrid(max(col)+1) - marg(1)/2;
w = xGrid(min(col)) + marg(1)/2;

% return position
pos = [w, s, abs(e-w), abs(n-s)];
end
%===========================================================