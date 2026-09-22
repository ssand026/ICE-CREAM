% Set up the ice-cream cone system
%-----------------------------------------------------------
clc; clearvars;

% load the mesh parameters
mesh = FE_mesh.load("ice_cream.step",15e3,"maxIters",20);
mesh = mesh.smooth("fixedSurfs","all","iters",10);

% define units
[nm,eV,fs] = atomic_units("nm","eV","fs");

% scale mesh
zScl = 1/(max(mesh.pts(3,:))-min(mesh.pts(3,:)));
mesh.pts = 100*nm * zScl * mesh.pts;
clearvars zScl

% re-order mesh
mesh = mesh.restructure("amd");

% display
meshplot3D(mesh.pts,mesh.tri,[],"meshStyle",{"FaceAlpha",1})


%% define system properties
%-----------------------------------------------------------
% set the inverse-mass tensor
m_xy = 0.19;
m_zz = 0.98;
massCoeff = diag(1./(2*[m_xy, m_xy, m_zz]));
massCoeff = FE_coeff(mesh,massCoeff);
clear m_xy m_zz

% set the spatially-dependent conduction-band edge
bandCoeff = "@(zz) max(zz)-zz";
bandCoeff = (50e-3*eV/nm) * FE_coeff(mesh,bandCoeff);

% assemble the Hamiltonian
%-----------------------------------------------------------
BCM = FEmat(mesh,"dirichlet");
M = FEmat(mesh,"overlap");

K = FEmat(mesh,"stiffness",massCoeff);
V = FEmat(mesh,"scalar",bandCoeff);
H = K + V;

[muE,muB,muZ] = FEmat_dipoles(mesh,massCoeff,-1);
clearvars K V massCoeff bandCoeff

% convert matrices to the reduced basis
H = to_basis(H,"reduced",BCM{:});
M = to_basis(M,"reduced",BCM{:});

muE = to_basis(muE,"reduced",BCM{:});
muB = to_basis(muB,"reduced",BCM{:});
muZ = to_basis(muZ,"reduced",BCM{:});

%% get eigenstates
%-----------------------------------------------------------
[states,E] = eigenstates(H,M,"numEigs",14);

% plot eigenstates
data = real(to_basis(states,"full",BCM{:}));
meshplot3D(mesh.pts,mesh.tri,data,"style","volume","numSurface",8, ...
	"cmapLimits",["-MAG/2","+MAG/2"],"faceReduce",0)

%% save results
%-----------------------------------------------------------
% save("ice_cream_sys.mat")










