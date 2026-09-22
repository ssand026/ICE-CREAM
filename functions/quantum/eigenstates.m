function [states, energy, didConvg] = eigenstates(H,M,opt)
% Returns the normalized eigenvectors and the eigenvalues for a Hamiltonian system.
% Specifying the overlap-matrix "M" as the 2nd argument solves for the eigen-
% values/vectors H*v = λ*M*v instead of H*v = λ*v
arguments
	H double {mustBeSquare}
	M double {mustBeSquare} = [];
	opt.numEigs  (1,1) double {mustBeInteger,mustBePositive} = 16;
	opt.tolerance (1,1) double {mustBeInRange(opt.tolerance,0,1)} = 8*eps;
	opt.maxIters  (1,1) double {mustBeInteger,mustBePositive} = 512;
	opt.isPosDef  (1,1) logical = false;
	opt.range (1,1) {mustBeMember(opt.range,["smallestreal","smallestabs",])} = "smallestreal";
end

% set eigensolver options
numEigs = min(opt.numEigs,length(H));
numSubspace = max(2*(numEigs+2),32);
opts = {numEigs, opt.range, "FailureTreatment","keep", "IsSymmetricDefinite",opt.isPosDef, ...
	"SubspaceDimension",numSubspace, "MaxIterations",opt.maxIters, "Tolerance",opt.tolerance};

if iseye(M)
	% orthonormal basis
	%-----------------------------
	% initialize start vector
	v0 = full(diag(H));
	v0 = v0 - min(v0);
	v0 = v0 + max(v0)/2;

	% get eigenstates
	[states,energy,flag] = eigs(H,opts{:},"StartVector",v0);
	states = states ./ sqrt(sum( abs(states).^2, 1));
else
	% non-orthonormal basis
	%-----------------------------
	% initialize start vector
	v0 = full(diag(M));
	v0 = v0 - min(v0);
	v0 = v0 + max(v0)/2;

	% get eigenstates
	[states,energy,flag] = eigs(H,M,opts{:},"StartVector",v0);
	states = states ./ sqrt(sum( abs(conj(states) .* (M * states)), 1));
end
energy = reshape(diag(energy),1,[]);
didConvg = ~flag;
% DONE
end