function [psi] = propagate_MH(psi, tau, H, M, iM, dM, opt)
% Forward-propagates a set of states using using the Al-Mohy/Higham exponentiation method.
%
% Their method for computing the action of a matrix exponential on a vector is described in:
% Awad H. Al-Mohy and Nicholas J. Higham. "Computing the action of the matrix exponential, 
% with an application to exponential integrators." SIAM J. Sci. Comput., 33(2):488-511, 2011.
% DOI. 10.1137/100788860
% 
% SEE ALSO: EXPM, EXPMV
arguments
	psi (:,:) double
	tau (1,1) double
	H  double {mustBeSquare}
	M  double {mustBeSquare} = [];
	iM double {mustBeSquare} = [];
	dM {mustBeA(dM,{'double','decomposition','preconditioned'})} = [];
	opt.minTerms (1,1) double {mustBeInteger,mustBePositive} = 3;
	opt.maxTerms (1,1) double {mustBeInteger,mustBePositive} = 150;
	opt.maxScale (1,1) double {mustBeInteger,mustBePositive} = 1000;
	opt.tolerance (1,1) double {mustBeInRange(opt.tolerance,0,1)} = 2*eps;
	opt.warnings (1,1) logical = false;
	opt.fastNorm (1,1) logical = false;
end

%% check/parse inputs
%===========================================================
% check that the inputs have appropriate sizes
if (size(H,2) ~= size(psi,1))
	error("ERROR: the sizes of 'H' and 'psi' must match")
elseif ~isempty(M) && ~isequal(size(M),size(H))
	error("ERROR: the sizes of 'H' and 'M' must match")
elseif ~isempty(iM) && ~isequal(size(iM),size(M))
	error("ERROR: the sizes of 'M' and 'iM' must match")
elseif opt.minTerms >= opt.maxTerms
	error("ERROR: the value of 'minTerms' must be less than 'maxTerms'")
elseif opt.minTerms > numel(theta_eps)
	error("ERROR: the value of 'minTerms' must be less than or equal to "+numel(theta_eps))
end

% check that the inputs are appropriately bound
opt.maxTerms = min(opt.maxTerms,numel(theta_eps));
opt.tolerance = max(opt.tolerance,eps);

% check for trivial input cases
if (abs(tau)==0)
	% exponent is zero-valued, return psi and exit
	return
elseif all(isnan(psi(:)) | (psi(:)==0))
	% psi will not change, return psi and exit
	return
end

%% apply shift to "H" for better Taylor-series convergence
%===========================================================
% the optimal shift "mu" should minimize the Frobenius norm of "H"

if iseye(M)
	% the basis is orthonormal
	%-----------------------------------------------------------
	if isdiag(H)
		% matrix is diagonal, return exponentiation and exit
		psi = exp(-1i*tau*diag(H)) .* psi;
		return
	end
	
	% set the decomposion of M
	dM = 1;

	% calculate shift from H directly
	lenH = length(H);
	mu = trace(H)/lenH;
	diagIndx = [1:(lenH+1):numel(H)];
	H(diagIndx) = H(diagIndx) - mu; % H = H - mu*I

	% get the 1-norm of H
	norm_H = norm(H,1);
else
	% the basis is not orthonormal
	%-----------------------------------------------------------
	% get the decomposion and inverse of M
	if isa(dM,"double")
		dM = decomposition(M);
	end
	if isempty(iM)
		iM = dM \ eye(size(M));
		iM = (iM + iM.')/2;
	end
	
	% compute shift from the Frobenius inner-product of inv(M) and H
	lenH = length(H);
	mu = (iM(:).' * H(:))/lenH; % mu = trace(M\H)/sqrt(numel(H))
	H = H - mu*M;  % H = H - mu*M

	% get the 1-norm of M\H
	if opt.fastNorm
		% estimate using a subset of the columns of H
		numCols = min(max(16,ceil(sqrt(lenH)/2)),lenH);
		[~,maxCols] = maxk(abs(diag(iM).*diag(H)),numCols);
		norm_H = norm(iM*H(:,maxCols),1);
	else
		% compute the full norm
		norm_H = norm(iM * H,1);
	end
end

%% deterine the optimal values for the scaling factor 's' and Taylor degree 'm'
%===========================================================
[s,m,cost] = scaling_params(abs(tau)*norm_H, opt.tolerance, opt.minTerms, opt.maxTerms);
maxCost = opt.maxScale * opt.maxTerms;
if (cost > maxCost) && opt.warnings
	% warn if exponentiation cost is large
	warning("The expected cost of the calculation ("+cost+") exceeds the given "+ ...
		"cost threshold ("+maxCost+"). Relax the solution tolerance or decrease "+ ...
		"the step-size to avoid this problem.")
end

%% compute the scaled taylor-series
%===========================================================
alpha = exp(-1i*tau*mu/s);
beta = -1i*tau/s;
for ss = 1:s
	% check if any columns of psi are nilpotent/non-finite
	converged = all(psi==0 | ~isfinite(psi),1);
	if all(converged); break; end

	% compute the taylor series, updating only the non-converged vectors
	t_m = psi;
	for mm = 1:m
		% get the norm of the previous update vector
		prevNorm = vecnorm(t_m,2,1);

		% obtain the next term in the taylor series
		t_m = (beta/mm) * (dM \ (H * t_m));
		psi = psi + t_m;

		% check if any vectors have converged
		currNorm = vecnorm(t_m,2,1);
		fullNorm = vecnorm(psi,2,1);
		converged = converged | ((prevNorm + currNorm) < opt.tolerance * fullNorm);

		% exit taylor-series expansion if all vectors have converged
		if all(converged); break; end
	end
	% apply scaling factor
	psi = alpha * psi;
end

% % the alternate version below acts on the non-converged vectors only. 
% % In theory it should require less computation, but indexing into the
% % vectors might prevent certain memory optimizations.
% %-----------------------------------------------------------
% for ss = 1:s
% 	% check if any columns of psi are nilpotent/non-finite
% 	converged = all(psi==0 | ~isfinite(psi),1);
% 	if all(converged); break; end
% 
% 	% compute the taylor series, updating only the non-converged vectors
% 	t_m = psi;
% 	for mm = 1:m
% 		% get the norm of the previous update vector
% 		prevNorm = vecnorm(t_m(:,~converged),2,1);
% 
% 		% obtain the next term in the taylor series
% 		t_m(:,~converged) = (beta/mm) * (dM \ (H * t_m(:,~converged)));
% 		psi(:,~converged) = psi(:,~converged) + t_m(:,~converged);
% 
% 		% check if any vectors have converged
% 		currNorm = vecnorm(t_m(:,~converged),2,1);
% 		fullNorm = vecnorm(psi(:,~converged),2,1);
% 		converged(~converged) = ((prevNorm + currNorm) < opt.tolerance * fullNorm);
% 
% 		% exit taylor-series expansion if all vectors have converged
% 		if all(converged); break; end
% 	end
% 	% apply scaling factor
% 	psi = alpha * psi;
% end

%% DONE
end

% external functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,m,cost] = scaling_params(normH, tol, minDegree, maxDegree)
% deterine the optimal values for the scaling factor 's' and taylor degree 'm'
m_vals = [minDegree : maxDegree];
theta_m = reshape(theta_eps,1,[]);
theta_m = theta_m(m_vals);

% scale by the m-th root of (tol/eps) to obtain theta_m for the specified tolerance
if (tol ~= eps)
	mroot = 2.^((log2(tol)-log2(eps))./m_vals); % mroot = nthroot(tol/eps, m_vals)
	theta_m = mroot .* theta_m;
end

% find the minimum value of 's' such that normH/s <= norm_thresh
m_max = max(m_vals);
p_max = floor(max(roots([1,-1,-m_max])));
norm_thresh = 4 * theta_m(end) * p_max * (p_max+3) / m_vals(end);
s_min = (1 + floor(normH/norm_thresh));

% find which values of 's' and 'm' minimize the computational cost 
s_vals = max(s_min, ceil(normH ./ theta_m));
[cost,index] = min(s_vals .* m_vals);
s = s_vals(index);
m = m_vals(index);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta] = theta_eps(~)
% first 240 values of theta_m given a tolerance of 2^(-52), i.e. the machine
% tolerance for double-typed floating point data
%
% theta_m was computed using the Mathematica (14.1) code:
%	pMax = 16;
%   mMax = pMax*(pMax - 1);
%   eps = 2^-52;
%
%   T[m_, x_] := Normal[Series[Exp[x], {x, 0, m}]]
%	c[m_, n_] := CoefficientList[Series[Log[T[m, x]] - x, {x, 0, n}], x];
%	h[m_, n_, x_] := Normal[SeriesData[x, 0, Abs[c[m, n]], -1, n - 1, 1]];
%
%	theta[m_, tol_] := 
%	  For[nn=(2*m+1); yOld=0; error=1;, (error >= 2^-52), nn++,
%		hN = h[m, nn, x];
%		yNew = N[Max[x /. Assuming[x \[Element] Reals, Solve[hN == tol, x]]]];
%		error = Abs[yNew - yOld]/(yNew + yOld);
%		yOld = yNew;
%		If[(error < 2^-52), Return[yNew]];
%	  ];
%
%   thetaList = Map[theta[#, eps] &, Range[mMax]]
%
% SEE ALSO: EPS, SCALING_PARAMS
theta = [...
	4.440892098500624e-16; 3.650024100028822e-08; 1.746687180598474e-05;
	4.039883260957763e-04; 2.757714014578019e-03; 1.017421010361573e-02;
	2.631825279562541e-02; 5.440223724700706e-02; 9.667973113502470e-02;
	1.543868760221968e-01; 2.279079904943047e-01; 3.170045301813742e-01;
	4.210278654598218e-01; 5.390880427645200e-01; 6.701759949669360e-01;
	8.132464267738090e-01; 9.672708106351640e-01; 1.131269058497103e+00;
	1.304326648308250e+00; 1.485602195818629e+00; 1.674328983539579e+00;
	1.869812840386788e+00; 2.071427961642975e+00; 2.278611699624979e+00;
	2.490858975313664e+00; 2.707716707539473e+00; 2.928778489814883e+00;
	3.153679637466287e+00; 3.382092659772554e+00; 3.613723169977843e+00;
	3.848306221266829e+00; 4.085603043018848e+00; 4.325398144876080e+00;
	4.567496753723782e+00; 4.811722548859454e+00; 5.057915662305048e+00;
	5.305930913678703e+00; 5.555636251842683e+00; 5.806911378409393e+00;
	6.059646530957151e+00; 6.313741406394577e+00; 6.569104207274844e+00;
	6.825650795985463e+00; 7.083303943629002e+00; 7.341992662078180e+00;
	7.601651609153164e+00; 7.862220558149172e+00; 8.123643924059080e+00;
	8.385870339807280e+00; 8.648852276656010e+00; 8.912545703678880e+00;
	9.176909781834280e+00; 9.441906588724410e+00; 9.707500870607390e+00;
	9.973659818648020e+00; 1.024035286675660e+01; 1.050755150868279e+01;
	1.077522913230716e+01; 1.104336086931504e+01; 1.131192345864726e+01;
	1.158089512230750e+01; 1.185025545226655e+01; 1.211998530734631e+01;
	1.239006671908980e+01; 1.266048280573346e+01; 1.293121769349372e+01;
	1.320225644446506e+01; 1.347358499050138e+01; 1.374519007251879e+01;
	1.401705918471649e+01; 1.428918052326429e+01; 1.456154293905141e+01;
	1.483413589413206e+01; 1.510694942153965e+01; 1.537997408817407e+01;
	1.565320096049517e+01; 1.592662157278145e+01; 1.620022789773614e+01;
	1.647401231924326e+01; 1.674796760709501e+01; 1.702208689352820e+01;
	1.729636365142253e+01; 1.757079167402660e+01; 1.784536505609011e+01;
	1.812007817629098e+01; 1.839492568085631e+01; 1.866990246828482e+01;
	1.894500367508637e+01; 1.922022466246145e+01; 1.949556100384994e+01;
	1.977100847328467e+01; 2.004656303449031e+01; 2.032222083067332e+01;
	2.059797817495290e+01; 2.087383154138722e+01; 2.114977755655250e+01;
	2.142581299163616e+01; 2.170193475500838e+01; 2.197813988523873e+01;
	2.225442554452773e+01; 2.253078901252513e+01; 2.280722768050876e+01;
	2.308373904590017e+01; 2.336032070709469e+01; 2.363697035858533e+01;
	2.391368578636150e+01; 2.419046486356483e+01; 2.446730554638570e+01;
	2.474420587018523e+01; 2.502116394582858e+01; 2.529817795621648e+01;
	2.557524615300256e+01; 2.585236685348542e+01; 2.612953843766446e+01;
	2.640675934544991e+01; 2.668402807401766e+01; 2.696134317530027e+01;
	2.723870325360644e+01; 2.751610696336095e+01; 2.779355300695863e+01;
	2.807104013272528e+01; 2.834856713297993e+01; 2.862613284219224e+01;
	2.890373613523005e+01; 2.918137592569183e+01; 2.945905116431944e+01;
	2.973676083748669e+01; 3.001450396575967e+01; 3.029227960252486e+01;
	3.057008683268141e+01; 3.084792477139416e+01; 3.112579256290422e+01;
	3.140368937939393e+01; 3.168161441990358e+01; 3.195956690929695e+01;
	3.223754609727338e+01; 3.251555125742378e+01; 3.279358168632853e+01;
	3.307163670269497e+01; 3.334971564653265e+01; 3.362781787836430e+01;
	3.390594277847095e+01; 3.418408974616927e+01; 3.446225819911972e+01;
	3.474044757266406e+01; 3.501865731919059e+01; 3.529688690752593e+01;
	3.557513582235210e+01; 3.585340356364767e+01; 3.613168964615168e+01;
	3.640999359884961e+01; 3.668831496448004e+01; 3.696665329906120e+01;
	3.724500817143650e+01; 3.752337916283810e+01; 3.780176586646784e+01;
	3.808016788709445e+01; 3.835858484066675e+01; 3.863701635394172e+01;
	3.891546206412699e+01; 3.919392161853705e+01; 3.947239467426265e+01;
	3.975088089785266e+01; 4.002937996500809e+01; 4.030789156028745e+01;
	4.058641537682323e+01; 4.086495111604884e+01; 4.114349848743570e+01;
	4.142205720823997e+01; 4.170062700325849e+01; 4.197920760459367e+01;
	4.225779875142686e+01; 4.253640018979977e+01; 4.281501167240394e+01;
	4.309363295837740e+01; 4.337226381310880e+01; 4.365090400804829e+01;
	4.392955332052510e+01; 4.420821153357138e+01; 4.448687843575224e+01;
	4.476555382100168e+01; 4.504423748846404e+01; 4.532292924234092e+01;
	4.560162889174340e+01; 4.588033625054909e+01; 4.615905113726409e+01;
	4.643777337488960e+01; 4.671650279079285e+01; 4.699523921658249e+01;
	4.727398248798793e+01; 4.755273244474275e+01; 4.783148893047187e+01;
	4.811025179258237e+01; 4.838902088215798e+01; 4.866779605385679e+01;
	4.894657716581230e+01; 4.922536407953777e+01; 4.950415665983343e+01;
	4.978295477469667e+01; 5.006175829523513e+01; 5.034056709558243e+01;
	5.061938105281666e+01; 5.089820004688120e+01; 5.117702396050813e+01;
	5.145585267914406e+01; 5.173468609087804e+01; 5.201352408637185e+01;
	5.229236655879229e+01; 5.257121340374562e+01; 5.285006451921390e+01;
	5.312891980549325e+01; 5.340777916513399e+01; 5.368664250288251e+01;
	5.396550972562498e+01; 5.424438074233255e+01; 5.452325546400823e+01;
	5.480213380363546e+01; 5.508101567612790e+01; 5.535990099828096e+01;
	5.563878968872457e+01; 5.591768166787730e+01; 5.619657685790188e+01;
	5.647547518266198e+01; 5.675437656768015e+01; 5.703328094009698e+01;
	5.731218822863142e+01; 5.759109836354233e+01; 5.787001127659080e+01;
	5.814892690100387e+01; 5.842784517143905e+01; 5.870676602394989e+01;
	5.898568939595242e+01; 5.926461522619265e+01; 5.954354345471487e+01;
	5.982247402283081e+01; 6.010140687308959e+01; 6.038034194924870e+01;
	6.065927919624546e+01; 6.093821856016943e+01; 6.121715998823562e+01];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%