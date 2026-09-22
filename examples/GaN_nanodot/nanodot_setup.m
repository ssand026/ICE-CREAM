%% Set up the GaN wurtztite nanodot
%-----------------------------------------------------------
clc; clearvars;

% define units
[nm,eV,fs] = atomic_units("nm","eV","fs");

% generate the initial mesh
mesh = FE_mesh.load("GaN_nanodot.step",11e3,"maxIters",20);

% smooth/restructure/rescale
mesh = mesh.restructure("amd");
mesh = mesh.smooth("iters",5,"fixedSurfs","all");
yScl = 1/(max(mesh.pts(2,:))-min(mesh.pts(2,:)));
mesh.pts = 12.5*nm * yScl * mesh.pts;
clearvars yScl

% find feature edges
faces = double(mesh.boundaryFaces);
faces = align_normals(faces);
warning("off","MATLAB:triangulation:PtsNotInTriWarnId")
TR = triangulation(faces,mesh.pts.');
warning("on","MATLAB:triangulation:PtsNotInTriWarnId")
edges = TR.featureEdges(15*pi/180);

xEdg = reshape(TR.Points(edges.',1),2,[]);
yEdg = reshape(TR.Points(edges.',2),2,[]);
zEdg = reshape(TR.Points(edges.',3),2,[]);
clearvars faces TR edges

% display
meshplot3D(mesh.pts,mesh.tri,[],"meshStyle",{"FaceAlpha",0.75,"EdgeAlpha",0.25})
hold on
plot3(xEdg,yEdg,zEdg,"-k") % highlight feature edges
hold off

%% define the system Hamiltonian
%-----------------------------------------------------------
% set the inverse-mass tensor
m_xy = 0.2;
m_zz = 0.2;
massCoeff = diag(1./(2*[m_xy, m_xy, m_zz]));
massCoeff = FE_coeff(mesh,massCoeff);
clear m_xy m_zz

% set the spatially-dependent conduction-band edge
z = mesh.pts(3,:);
zNorm = (z - min(z))/(max(z)-min(z));
bandCoeff = (4 - zNorm) * eV;
bandCoeff = FE_coeff(mesh,bandCoeff);
clearvars z zNorm

% assemble components of the Hamiltonian
BCM = FEmat(mesh,"dirichlet");
M = FEmat(mesh,"overlap");
K = FEmat(mesh,"stiffness",massCoeff);
V = FEmat(mesh,"scalar",bandCoeff);
[muE,muB,muZ] = FEmat_dipoles(mesh,massCoeff);

H = K + V; % Hamiltonian
clearvars K V

% convert matrices to the reduced basis
H = to_basis(H,"reduced",BCM{:});
M = to_basis(M,"reduced",BCM{:});
muE = to_basis(muE,"reduced",BCM{:});
muB = to_basis(muB,"reduced",BCM{:});
muZ = to_basis(muZ,"reduced",BCM{:});

%% get eigenstates
%-----------------------------------------------------------
[states,E] = eigenstates(H,M,"numEigs",16);

% plot eigenstates
meshplot3D(mesh.pts,mesh.tri,to_basis(states,"full",BCM{:}),"numSurface",12, ...
	"style","surface","cmapLimits",["-2.5*STD","+2.5*STD"],"axesStyle","arrows")
hold on
plot3(xEdg,yEdg,zEdg,"-","Color",0*[1 1 1])
hold off


% check coupling between states via the electric/magnetic fields
%-----------------------------------------------------------
cmap = cmap_mono;
cmap = cmap(3:end,:);
cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,40).^2);

TDM_E = vecnorm(get_TDM(states,muE(1:2),"noDiags",true),2,3);
TDM_E = log10(TDM_E);
figure; imagesc(TDM_E); axis image
clim(max(TDM_E(:))+[-4,0])
colorbar;
colormap(cmap)

TDM_B = vecnorm(get_TDM(states,muB(1:2),"noDiags",true),2,3);
TDM_B = log10(TDM_B);
figure; imagesc(TDM_B); axis image
clim(max(TDM_B(:))+[-4,0])
colorbar;
colormap(cmap)

%% save results
%-----------------------------------------------------------
save("nanodot_sys.mat")
