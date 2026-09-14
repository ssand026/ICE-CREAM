%% load the nanodot system parameters
clc; clearvars;
load("nanodot_sys.mat")

% run control using electric fields
mu = muE;

% determine starting/target states/occupancies
ni = 4;
nf = 7;
psi_i = states(:,ni);
psi_f = states(:,nf);

% get transition diagnostics
[muRot,rotT] = rotate_dipoles(psi_i,psi_f,mu);
QSL = get_QSL("fermi",states,H,M);
TDM = vecnorm(get_TDM(states,mu),2,3);

% determine the timestep/span values
numQSL = 20;		% multiply by the QSL to obtain the total time
tMult = 10 * fs;	% ensure total duration is a multiple of tMult

T = ceil(numQSL * QSL(ni,nf)/tMult)*tMult;
numTau = 200 * numQSL;
tgrid = linspace(0,T,numTau+1);

fprintf("T = %i fs \n",round(T/fs));
clearvars numQSL tMult

%% test the pulse-area-theorem (PAT) pulse
%-----------------------------------------------------------

% compute the PAT-predicted pulse
env = pulse_envelope(tgrid,"sine",[T/2,T/2]);
envArea = (env(1:end-1)+env(2:end))/2 * diff(tgrid(:));

fldDir = rotT(:,1);
fldMag = pi./(envArea * TDM(ni,nf));
fldFreq = 1/QSL(ni,nf);
field = fldDir .* fldMag .* (env .* sin(2*pi*fldFreq*tgrid));
clearvars envArea

% set optimization algorithm options
QC_opts = { ...
	"linesearch",	0, ...
	"maxIters",		0, ...
	"plotting",		true, ...
	"propagator",	"CN", ...
	"propOpts",		{"method","cgs"}, ...
	"envelope",		1, ...
	"stallTol",		0, ...
	};

[P,Pt,dP,~,psiFwd] = optimize_field(psi_i,psi_f',tgrid,field,H,mu,M,QC_opts{:});

% save
save("nanodot_4_7_PAT.mat")