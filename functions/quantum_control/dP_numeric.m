function [dP] = dP_numeric(psiFwd, phiBwd, tgrid, field, H, mu, M, iM, dM, opt)
% DP_NUMERIC computes the change in the state-to-state transition probability with respect to a
% time-dependent control-field with dipole moments "mu" using a central differences approach
arguments
	psiFwd (:,:,:) double
	phiBwd (:,:,:) double
	tgrid (1,:) double
	field (:,:) double
	H  double {mustBeSquare}
	mu (:,:) cell
	M  double {mustBeSquare} = [];
	iM double {mustBeSquare} = [];
	dM {mustBeA(dM,{'double','decomposition'})} = [];
	% propagator options
	opt.propagator (1,1) {mustBeMember(opt.propagator,["MH","CN"])}
	opt.propOpts (1,:) cell = {};
	opt.normalize (1,1) logical = true;
	% options for self-consistent field calculations
	opt.basisM (1,2) cell
	opt.PoissonEq PoissonSolver
	opt.occ_i (:,1) {mustBeFinite,mustBeNonnegative} = 1;
	opt.occ_f (:,1) {mustBeFinite,mustBeNonnegative} = 1;
	opt.minScfIters (1,1) {mustBeInteger,mustBeNonnegative} = 2;
	opt.maxScfIters (1,1) {mustBeInteger,mustBeNonnegative} = 12;
	opt.scfTol (1,1) {mustBePositive,mustBeInRange(opt.scfTol,0,1)} = 1e-6;
end

% configure the propagator
%-----------------------------------------------------------
% check if the basis is orthonormal
isOrth = iseye(M);

% intialize the decomposed/inverted overlap matrix if needed
if ~isOrth && (opt.propagator=="MH")
	if isa(dM,"double")
		dM = decomposition(M);
	end
	if isempty(iM)
		iM = dM \ eye(size(M));
		iM = (iM + iM.')/2;
	end
else
	dM = [];
	iM = [];
end

% set basic options
propArgs = {H, mu, M, iM, dM,"propagator",opt.propagator,"propOpts",opt.propOpts, ...
	"normalize",opt.normalize,"display",false,"saveAll",false};
fwdArgs = {"occupancy",opt.occ_i,"direction","fwd"};
bwdArgs = {"occupancy",opt.occ_f,"direction","bwd"};

% check if self-interaction is enabled
scfEnabled = isfield(opt,"PoissonEq") && isfield(opt,"basisM");
if scfEnabled
	% set SCF options
	propArgs = [propArgs,{"PoissonEq",opt.PoissonEq,"basisM",opt.basisM,"scfTol",opt.scfTol}];
	propArgs = [propArgs,{"minScfIters",opt.minScfIters,"maxScfIters",opt.maxScfIters}];
end

% compute the gradient terms
%-----------------------------------------------------------
dP = cell(1,numel(mu));
numPhi = size(phiBwd,1);
numPsi = size(psiFwd,2);
muIndex = find(~cellfun(@isempty,mu));

for dd = muIndex
	% determine a minimum shift "h" so that all updates to the Hamiltonian
	% by the dipole matrix "mu" fall within numerical precision
	nonZero = (mu{dd}~=0);
	mu_vals = mu{dd}(nonZero);
	H_vals = H(nonZero);
	h_min = max(eps(H_vals)./abs(mu_vals),[],"all");

	% compute the gradient values
	%-----------------------------------------------------------
	dP{dd} = zeros(numPhi,numPsi,numel(tgrid));
	for tt = 1:numel(tgrid)
		% determine "h" for the current field value
		field_tol = eps(field(dd,tt));
		h = 2 * max(h_min,field_tol);
		h_dd = zeros(size(field,1),1);
		h_dd(dd) = h;

		if tt == 1
			% use initial psi
			psi_pos = psiFwd(:,:,tt);
			psi_neg = psiFwd(:,:,tt);
		else
			% compute psi(:,:,tt) when field(dd,tt) is shifted by +h/-h
			psi_pos = propagate( psiFwd(:,:,tt-1), tgrid(:,tt-1:tt), ...
				[field(:,tt-1),field(:,tt)+h_dd], propArgs{:}, fwdArgs{:});
			psi_neg = propagate( psiFwd(:,:,tt-1), tgrid(:,tt-1:tt), ...
				[field(:,tt-1),field(:,tt)-h_dd], propArgs{:}, fwdArgs{:});
		end
		
		if tt==numel(tgrid)
			% use final phi
			phi_pos = phiBwd(:,:,tt);
			phi_neg = phiBwd(:,:,tt);
		else
			% compute phi(:,:,tt) when field(dd,tt) is shifted by +h/-h
			phi_pos = propagate( phiBwd(:,:,tt+1), tgrid(:,tt:tt+1), ...
				[field(:,tt)+h_dd,field(:,tt+1)], propArgs{:}, bwdArgs{:});
			phi_neg = propagate( phiBwd(:,:,tt+1), tgrid(:,tt:tt+1), ...
				[field(:,tt)-h_dd,field(:,tt+1)], propArgs{:}, bwdArgs{:});
		end

		% compute the gradient values
		if isOrth
			P_pos = abs(phi_pos * psi_pos).^2;
			P_neg = abs(phi_neg * psi_neg).^2;
		else
			P_pos = abs(phi_pos * M * psi_pos).^2;
			P_neg = abs(phi_neg * M * psi_neg).^2;
		end
		dP{dd}(:,:,tt) = (P_pos - P_neg)/(2*h);
	end
end
% DONE
end