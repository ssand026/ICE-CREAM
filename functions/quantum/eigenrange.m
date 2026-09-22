function [eigVec, eigVal] = eigenrange(limits,H,M,opt)
% Returns the eigenvectors and the eigenvalues within the given range of values
% Specifying the overlap-matrix "M" as the 3rd argument solves for the eigen-
% values/vectors H*v = λ*M*v instead of H*v = λ*v
arguments
	limits (1,2) double {mustBeReal}
	H double {mustBeSquare}
	M double {mustBeSquare} = [];
	opt.maxEigs (1,1) double {mustBeInteger,mustBePositive} = intmax;
	opt.tolerance (1,1) double {mustBeInRange(opt.tolerance,0,1)} = 1e-12;
end
if limits(1)>=limits(2)
	error("ERROR: the input 'limits' must contain ascending values")
end

% parse options
basisSize = size(H,1);
maxEigs = min(opt.maxEigs,basisSize);
relError = @(x) 2*eps(x)*(opt.tolerance/eps(1));
eigsOptions = {"FailureTreatment","keep","Tolerance",opt.tolerance/2};

% check if the overlap matrix is the identity matrix
[isEye,scale] = iseye(M);
if ~isEye && isfinite(scale)
	% the overlap can be converted to an identity matrix
	H = H / scale;
	isEye = true;
end

% set up the iterative eigensolver
if isEye
	localEig = @(sigma,numEig) eigs(H, numEig, sigma, ...
		"SubspaceDimension", 2*ceil(1.5*(numEig+1)), eigsOptions{:});
else
	localEig = @(sigma,numEig) eigs(H, M, numEig, sigma, ...
		"SubspaceDimension", 2*ceil(1.5*(numEig+1)), eigsOptions{:});
end

% find the smallest eigenvalue given a non-finite lower bound
if (limits(1) == -Inf)
	[~,limits(1),~] = localEig("smallestreal",1);
end
% find the largest eigenvalue given a non-finite upper bound
if (limits(2) == +Inf)
	[~,limits(2),~] = localEig("largestreal",1);
end
% slightly expand limit window 
limits = limits + [-1,+1] .* relError(limits);


% begin extracting eigenvalues
%-----------------------------------------------------------
EIGS_PER_ITER = 6;
eigVal = zeros(1,0);
eigVec = zeros(basisSize,0);
sigma = limits(1); % center of the eigenvalue search

maxIters = ceil(2*maxEigs/EIGS_PER_ITER);
for iter = 1:maxIters
	if sigma > limits(2) || numel(eigVal) >= maxEigs
		% all eigenvals have been found
		break
	else
		% find additional eigenvals
		[vecs,vals,~] = localEig(sigma,EIGS_PER_ITER);
		vals = diag(vals);
	end

	% check for distinct eigenvals
	if ~isempty(eigVal)
		test = all(abs(eigVal-vals)>relError(eigVal),2);
		valid = vals(test);
	else
		valid = vals;
	end

	% check that these eigenvals fall within the search limits
	valid = valid(valid>=limits(1));
	valid = valid(valid<=limits(2));
	
	% add new values to the list
	newEigs = find(ismember(vals,valid));
	numNew = numel(newEigs);
	if numNew~=0
		eigVal(:,end+1:end+numNew) = vals(newEigs);
		eigVec(:,end+1:end+numNew) = vecs(:,newEigs);
	end

	% update sigma
	avgDiff = (max(eigVal)-min(eigVal))/numel(eigVal);
	sigma = max(sigma,max(eigVal)) + avgDiff;
end

% sort the eigenvectors and valus
[eigVal,order] = sort(eigVal(:));
eigVec = eigVec(:,order);

numEigs = numel(eigVal);
eigVal = eigVal(1:min(numEigs,maxEigs));
eigVec = eigVec(:,1:min(numEigs,maxEigs));

% normalize the eigenvectors
eigVec = norm_wavefunc(eigVec,M);

%DONE
end