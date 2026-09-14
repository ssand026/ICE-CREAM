%% load the ice-cream system parameters
clc; clearvars;
load("ice_cream_sys.mat")

% run control using magnetic fields
mu = muB;

% check coupling between states via the magnetic fields
%-----------------------------------------------------------
TDM = vecnorm(get_TDM(states,mu,"noDiags",true),2,3);

% display the coupling array
% figure; imagesc(TDM); axis image

% find the first several strongest couplings
%-----------------------------------------------------------
num = 5; % number of couplings
upperIndex = find(triu(ones(size(TDM))));
[strength,maxIndex] = maxk(TDM(upperIndex),num);
[ni,nf] = ind2sub(size(TDM),upperIndex(maxIndex)); % get the starting/ending state combinations
strength = strength/max(strength); % convert to relative coupling strengths
couplingList = [ni,nf,strength];
clearvars num ni nf upperIndex maxIndex strength

% disp(couplingList)

%% run QC for states with strong magnetic coupling
%-----------------------------------------------------------

% determine starting/target states/occupancies
psi_i = states(:,1);
psi_f = states(:,[6,7])*[1;-1i]/sqrt(2);

% set the timestep/simulation duration
QSL = get_QSL("fermi",[psi_i,psi_f],H,M);
T = 14 * QSL(1,2);
T = ceil(T/(1000*fs))*(1000*fs); % round to multiple of 1000 fs
fprintf("T = %i fs \n",round(T/fs));


numTau = 5000;
tgrid = linspace(0,T,numTau+1);
field = zeros(numel(mu),numTau+1);

% set optimization algorithm options
QC_opts = { ...
	"linesearch",	true, ...
	"minIters",		2, ...
	"maxIters",		3, ...
	"cutoff",		0.975, ...
	"plotting",		true, ...
	"propagator",	"MH", ...
	"propOpts",		{"fastNorm",true}, ...
	"envelope",		pulse_envelope(tgrid,"sine",[T/6,T/6]), ...
	"savefile",		"ice_cream_QC_prog.mat", ...
	};

% run the multi-state quantum control algorithm
[P,Pt,dP,field,psiFwd,exitFlag] = optimize_field(psi_i,psi_f',tgrid,field,H,mu,M,QC_opts{:});

% save("ice_cream_QC.mat")