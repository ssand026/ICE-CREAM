function [psi_t] = propagate(psi, tgrid, field, H, mu, M, iM, dM, opt)
% Computes the time-evolution of a set of states 'psi' according to the
% time-dependent Hamiltonian H(t) = H + field(t)*mu
%
% SEE ALSO: PROPAGATE_CN, PROPAGATE_MH, POISSONSOLVER
arguments (Input)
	psi   (:,:) double {mustBeFinite}
	tgrid (1,:) double {mustBeFinite}
	field (:,:) double {mustBeFinite}
	H  double {mustBeSquare}
	mu (1,:) cell
	M  double {mustBeSquare} = [];
	iM double {mustBeSquare} = [];
	dM {mustBeA(dM,{'double','decomposition','preconditioned'})} = [];
	% general options for propagate
	opt.propagator (1,1) {mustBeMember(opt.propagator,["MH","CN"])}
	opt.direction  (1,1) {mustBeMember(opt.direction,["fwd","bwd"])} = "fwd";
	opt.propOpts   (1,:) cell = {};
	opt.normalize (1,1) logical = true;
	opt.saveAll (1,1) logical = true;
	opt.display (1,1) logical = true;
	opt.progBar (1,1) {mustBeMember(opt.progBar,["simple","fancy","explicit"])} = "explicit";
end
arguments (Output)
	psi_t (:,:,:) double
end

%% prepare/check inputs
%===========================================================
% restructure the inputs if backwards-propagating
if (opt.direction == "bwd")
	tgrid = flip(tgrid,2);
	field = flip(field,2);
	psi = psi';
end

% check if inputs have valid size
tau = diff(tgrid,1,2);
numTau = numel(tau);
[numRow, numCol] = size(psi);
if (size(field,1) ~=  numel(mu))
	error("ERROR: the number of field components does not match the length of 'mu'.")
elseif ~isempty(field) && (size(field,2) ~= numel(tgrid))
	error("ERROR: the number of timesteps does not equal the number of field values.")
elseif (size(H,1) ~= numRow)
	error("ERROR: the sizes of 'H' and 'psi' must agree.")
elseif ~isempty(M) && ~isequal(size(M),size(H))
	error("ERROR: the sizes of 'H' and 'M' must match.")
elseif ~isempty(iM) && ~isequal(size(iM),size(M))
	error("ERROR: the sizes of 'M' and 'iM' must match.")
end

%% setup propagator
%===========================================================
% check if the basis is orthonormal
if iseye(M)
	propMethod = opt.propagator + "_orth";
	M = [];
else
	propMethod = opt.propagator;
end

% set propagation function
propOpt = opt.propOpts;
switch propMethod
	case "CN"
		% Crank-Nicolson propagation
		switch opt.direction
			case "fwd"; propFunc = @(x,t,HH) propagate_CN(x,t,HH,M,propOpt{:});
			case "bwd"; propFunc = @(x,t,HH) propagate_CN(x,t,HH.',M,propOpt{:});
		end
	case "CN_orth"
		% Crank-Nicolson propagation (orthonormal)
		I = eye(size(H),"like",H);
		switch opt.direction
			case "fwd"; propFunc = @(x,t,HH) propagate_CN(x,t,HH,I,propOpt{:});
			case "bwd"; propFunc = @(x,t,HH) propagate_CN(x,t,HH.',I,propOpt{:});
		end
	case "MH"
		% Al-Mohy/Higham propagation
		if isa(dM,"double")
			dM = decomposition(M);
		end
		if isempty(iM)
			iM = dM \ eye(size(M));
			iM = (iM + iM.')/2;
		end
		switch opt.direction
			case "fwd"; propFunc = @(x,t,HH) propagate_MH(x,t,HH,M,iM,dM,propOpt{:});
			case "bwd"; propFunc = @(x,t,HH) propagate_MH(x,t,HH.',M,iM,dM,propOpt{:});
		end
	case "MH_orth"
		% Al-Mohy/Higham propagation (orthonormal)
		switch opt.direction
			case "fwd"; propFunc = @(x,t,HH) propagate_MH(x,t,HH,propOpt{:});
			case "bwd"; propFunc = @(x,t,HH) propagate_MH(x,t,HH.',propOpt{:});
		end
end


%% pre loop processes
%===========================================================
if opt.normalize
	[psi,orig_scl] = norm_wavefunc(psi,M);
end
if opt.saveAll
	psi_t = complex(zeros(numRow,numCol,numTau+1));
	psi_t(:,:,1) = psi;
end
if opt.display
	if (opt.progBar == "explicit")
		dispBar = progbar("progress: ",opt.progBar,numTau);
	else
		dispBar = progbar("progress: ",opt.progBar,40);
	end
	dispBar(0);
end

%% propagation loop
%===========================================================
field_mid = (field(:,1:end-1) + field(:,2:end))/2;
for tt = 1:numTau
	% add time-dependent components to the Hamiltonian
	Ht = H;
	for dd = 1:numel(mu)
		Ht = Ht + field_mid(dd,tt) * mu{dd};
	end

	% apply propagator
	psi = propFunc(psi,tau(tt),Ht);
	if opt.normalize
		psi = norm_wavefunc(psi,M);
	end

	% add propagated wavefunctions to the output
	if opt.saveAll
		psi_t(:,:,tt+1) = psi;
	end
	
	% report propagation progress
	if opt.display
		dispBar(tt/numTau);
	end
end

%% post loop processes
%===========================================================
if (opt.saveAll == false)
	% return only the fully-evolved state
	psi_t = psi;
end
if opt.normalize
	% restore the scaling factor
	psi_t = psi_t .* orig_scl;
end

% restructure the output if performing backwards-propagation
if (opt.direction == "bwd")
	psi_t = conj(permute(psi_t,[2,1,3]));
	psi_t = flip(psi_t,3);
end
%% DONE
end