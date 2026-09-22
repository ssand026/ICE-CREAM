function [dP] = dP_analytic(psiFwd, phiBwd, tgrid, mu, M, opt)
% DP_ANALYTIC computes the change in the state-to-state transition probability with respect to a
% time-dependent control-field with dipole moments "mu".
arguments
	psiFwd (:,1,:) double
	phiBwd (1,:,:) double
	tgrid (1,:) double
	mu (:,:) cell
	M  double {mustBeSquare} = [];
	opt.propagator (1,1) {mustBeMember(opt.propagator,["MH","CN"])}
	opt.rephase (1,1) logical = false;
end

% set coefficients for gradient terms
switch opt.propagator
	case "CN"
		coeff = [1/4, 1/4, 1/4, 1/4];
	case "MH"
		coeff = [1/3, 1/6, 1/6, 1/3];
end

% get timestep duration/length
tau = reshape(diff(tgrid),1,1,[]);
numTau = numel(tau);

% get the expecation value for the overlap
if iseye(M)
	% basis is orthonormal
	P_total = conj(phiBwd(:,:,end) * psiFwd(:,:,end));
else
	% basis is non-orthonormal
	P_total = conj(phiBwd(:,:,end) * M * psiFwd(:,:,end));	
end


% compute the expectation values of mu at different propagation points
tauPos = ((1:numTau)+1); % index for the tau_{n+1/2} terms
tauNeg = ((1:numTau)+0); % index for the tau_{n-1/2} terms

dP = zeros(numel(mu),numel(tgrid));
muIndex = find(~cellfun(@isempty,mu));

for dd = muIndex
	% compute the gradient of the respective dipole moment
	mu_psi = tensormult(mu{dd},'ik',psiFwd,'kjt','ijt');

	dP_equal = pagemtimes(phiBwd, mu_psi); % 〈phi_{n}|mu|psi_{n}〉
	dP_a = dP_equal(:,:,tauPos);                               % 〈phi_{n+1}|mu|psi_{n+1}〉
	dP_b = pagemtimes(phiBwd(:,:,tauPos), mu_psi(:,:,tauNeg)); % 〈phi_{n+1}|mu|psi_{ n }〉
	dP_c = pagemtimes(phiBwd(:,:,tauNeg), mu_psi(:,:,tauPos)); % 〈phi_{ n }|mu|psi_{n+1}〉
	dP_d = dP_equal(:,:,tauNeg);                               % 〈phi_{n-1}|mu|psi_{n-1}〉

	% compute dP
	dP_mu = tau .* (coeff(1)*dP_a + coeff(2)*dP_b + coeff(3)*dP_c + coeff(4)*dP_d);
	dP_mu = cat(3,dP_mu(:,:,1), dP_mu(:,:,1:numTau-1)+dP_mu(:,:,2:numTau), dP_mu(:,:,numTau));
	dP(dd,:) = imag(P_total .* dP_mu);
end
% DONE
end