function [P,Pt,dP,field,psiFwd,exitFlag] = optimize_field(psi,phi,tgrid,field,H,mu,M,opt)
% OPTIMIZE_FIELD Optimize a time-dependent control field using gradient ascent.
%
% [Pt,dP,field,psiFwd,exitFlag] = OPTIMIZE_FIELD(psi,phi,tgrid,field,H,mu,M,opt) 
% performs gradient-based optimization of one or more control fields to maximize 
% the overlap between propagated initial states and a set of target states.
%
% The optimization repeatedly propagates the initial states forward in time, 
% computes the objective function, evaluates the gradient using either an 
% analytic adjoint method or finite differences, optionally performs a line 
% search to determine the optimal update step, and updates the control field 
% until convergence.
%
% The time-dependent Hamiltonian has the form:
%		H(t) = H - Σ_k μ_k E_k(t)
% where H is the field-free Hamiltonian, μ_k is the kth dipole-moment 
% operator, and E_k(t) is the kth-component of the control-field.
%
% Input Arguments
% ================
% psi   - Matrix of initial state vectors.
% phi   - Matrix of target state vectors.
% tgrid - Time grid used for propagation.
% field - Initial control field. Each row corresponds to a dipole operator in 'mu'.
% H     - Field-free Hamiltonian matrix.
% mu    - Cell array of dipole-moment operators.
% M     - Basis overlap matrix. Use EYE(N) or [] for an orthonormal basis.
% 
% Name-Value Arguments
% =====================
% propagator - propagation method ("MH" or "CN").
% learnRate  - Step-size when linesearch is disabled 
% linesearch - Enable adaptive line search.
% envelope   - fixed envelope for the field		 
% cutoff     - Exits if the overlap exceeds this value.
% minIters   - Minimum optimization iterations.
% maxIters   - Maximum optimization iterations.
% stallMax   - Exit if the overlap does not improve after this many iterations.
% stallTol   - Improvement threshold for the overlap.
%
% Output Arguments
% =================
% Pt       - Time-dependent overlap.
% dP       - Gradient of the overlap with respect to the control field
% field    - Optimized control field.
% psiFwd   - Forward-propagated wavefunctions corresponding to the optimized field.
% exitFlag - Optimization termination condition:
%               0: Objective threshold reached
%               1: Optimization process converged
%               2: Maximum iterations reached
%               
% 
% SEE ALSO: PROPAGATE, LINE_SEARCH, OVERLAP
arguments
	psi   (:,:) double
	phi   (:,:) double
	tgrid (1,:) double
	field (:,:) double
	H     (:,:) double {mustBeSquare}
	mu    (1,:) cell
	M     (:,:) double {mustBeSquare}
	% general options
	opt.propagator (1,1) {mustBeMember(opt.propagator,["MH","CN"])}
	opt.propOpts (1,:) cell = {};
	opt.savefile (1,:) string {mustBeScalarOrEmpty} = [];
	opt.display (1,1) logical = true;
	opt.plotting (1,1) logical = true;
	opt.parallel (1,1) logical = false;
	% optimization options
	opt.learnRate (1,1) double {mustBePositive} = 1;
	opt.linesearch (1,1) logical = false;
	opt.searchOpts (1,:) cell = {};
	opt.envelope (:,:) double {mustBeInRange(opt.envelope,0,1)} = 1;
	%
	opt.cutoff (1,1) double {mustBeInRange(opt.cutoff,0,1)} = 0.95;
	opt.minIters (1,1) uint16 {mustBeInteger,mustBeNonnegative} = 1;
	opt.maxIters (1,1) uint16 {mustBeInteger,mustBeNonnegative} = 10;
	opt.stallMax (1,1) uint16 {mustBeInteger,mustBePositive} = 2;
	opt.stallTol (1,1) {mustBeInRange(opt.stallTol,0,1)} = 0.005;
end

%% check inputs
%===========================================================
% transpose to form proper bra/kets
if (size(psi,1)~=length(H)) && (size(psi,2)==length(H)) 
	psi = psi';
end
if (size(phi,1)==length(H)) && (size(phi,2)~=length(H)) 
	phi = phi';
end

% check the size of the inputs
if (size(psi,1) ~= size(H,2))
	error("ERROR: the number of rows in 'psi' must equal the number of columns in 'H'.")
end
if (size(phi,2) ~= size(H,1))
	error("ERROR: the number of columns in 'phi' must equal the number of rows in 'H'.")
end
if ~iseye(M) && ~isequal(size(M),size(H))
	error("ERROR: the sizes of 'H' and 'M' must be equal.")
end
if ~isempty(field) && ~isequal(size(field),[numel(mu),numel(tgrid)])
	error("ERROR: the size of 'field' must equal be equal to [numel('mu') x numel('tgrid')].");
end
if ~all( cellfun(@(x)isequal(size(x),size(H)),mu,"UniformOutput",true) )
	error("ERROR: the sizes of 'mu' and 'H' must be equal.")
end
if all(size(opt.envelope,1)~=[size(field,1),1]) 
	error("ERROR: the pulse envelope must have rows equal to the number of field components.")
end
if all(size(opt.envelope,2)~=[numel(tgrid),1]) 
	error("ERROR: the pulse envelope must have columns equal to the number of timesteps.")
end

% check minIters
opt.minIters = min(opt.minIters,opt.maxIters);

% check parallelism
if (opt.parallel==true) && isempty(gcp("nocreate"))
	opt.parallel = false;
end

% check if saveProgress is enabled
if ~isempty(opt.savefile) && (strlength(opt.savefile) > 0)
	% check the savefile name
	[path,name,~] = fileparts(opt.savefile);
	if isvarname(name)
		saveProgress = true;
	else
		error("ERROR: invalid save-file name.")
	end
	
	% generate savefile path
	if endsWith(path,"/") || strlength(path)==0
		opt.savefile = path+name+".mat";
	else
		opt.savefile = path+"/"+name+".mat";
	end
	clearvars path name
else
	saveProgress = false;
end

%% configure the propagator
%===========================================================
% check if the basis is orthonormal
if ~iseye(M) && (opt.propagator=="MH")
	% intialize the decomposed/inverted overlap matrix
	dM = decomposition(M);
	% dM = preconditioned(M,"method","cgs");
	iM = dM \ eye(size(M));
	iM = (iM + iM.')/2;
else
	dM = [];
	iM = [];
end

propArgs = {H, mu, M, iM, dM,"propagator",opt.propagator,"propOpts",opt.propOpts,"normalize",true};

%% setup progress plots
%===========================================================
if opt.plotting
	fig = figure;
	fig.Position(3:4) = fig.Position(3) * [1, 1];
	tNorm = (0:numel(tgrid)-1);

	% field
	ax1 = axes("Parent",fig,"Units","normalized","NextPlot","ReplaceChildren", ...
		"PositionConstraint","outerposition","OuterPosition",[0, 2/3, 1, 1/3]);
	title(ax1,"\epsilon(t)","Interpreter","tex")
	xlim(tNorm([1,end]));
	% P(t)
	ax2 = axes("Parent",fig,"Units","normalized","NextPlot","ReplaceChildren", ...
		"PositionConstraint","outerposition","OuterPosition",[0, 1/3, 1, 1/3]);
	title(ax2,"P(t)","Interpreter","tex")
	xlim(tNorm([1,end]));
	% dP(t)
	ax3 = axes("Parent",fig,"Units","normalized","NextPlot","ReplaceChildren", ...
		"PositionConstraint","outerposition","OuterPosition",[0, 0/3, 1, 1/3]);
	title(ax3,"dP(t)","Interpreter","tex")
	xlim(tNorm([1,end]));
	%
	drawnow;
end


%% OPTIMIZATION LOOP
%===========================================================
% setup savefile
if saveProgress
	options = opt;
	save(opt.savefile,"tgrid","options","-v7.3");
	saveStruct = matfile(opt.savefile,"Writable",true);
end

exitFlag = -1;
for iter = 1:(1+opt.maxIters)
	%% LOOP START
	if opt.display
		fprintf("\n_____ Iteration %i _____\n",iter);
	end
	
	% plot the time-dependent control-field
	if opt.plotting && isgraphics(ax1)
		plot(tNorm,field.',"Parent",ax1);
		drawnow;
	end
	
	%% forwards propagate
	%===========================================================
	if opt.display; fprintf("forwards propagation: \n"); end
	%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	propCost = cputime;
	propTime = 0;
	tic;
	%...........................................................
	psiFwd = propagate(psi,tgrid,field,propArgs{:},"direction","fwd","display",opt.display,"saveAll",true);
	%...........................................................
	propTime = propTime + toc;
	propCost = cputime - propCost;
	%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	% save 
	if saveProgress
		saveStruct.field(1,iter) = {field};
		saveStruct.propTime(1,iter) = {propTime};
		saveStruct.propCost(1,iter) = {propCost};
	end

	%% compute the time-dependent overlap
	%===========================================================
	Pt = overlap(phi,M,psiFwd);
	P = Pt(:,:,end);

	% display progress
	if opt.display
		fprintf("P = %#4.5g \n",P);
	end
	
	% save progress
	if saveProgress
		saveStruct.P(1,iter) = {P};
		saveStruct.Pt(1,iter) = {Pt};
	end
	
	% plot P(t)
	if opt.plotting && isgraphics(ax2)
		Pt_max = max(Pt(:));
		Pt_max = log10(~Pt_max + Pt_max);
		Pt_max = 10^ceil(Pt_max);
		
		plot(tNorm, Pt(:).',"Parent",ax2);
		ax2.YLim = [0, Pt_max];
		clearvars Pt_max
		drawnow;
	end
	
	%% check exit criteria
	%===========================================================
	% evaluate the number of stagnant iterations
	if (iter == 1) || ((P - P_prev) > opt.stallTol)
		numStall = 0;
	else
		numStall = numStall + 1;
	end
	P_prev = P;

	% exit if ...
	if (iter <= opt.minIters)
		% do not exit
	elseif (opt.maxIters == 0)
		% exit after computing gradient
	elseif (P >= opt.cutoff)
		% the overlap meets the required threshold
		exitFlag = 0;
		break
	elseif (numStall >= opt.stallMax)
		% the field optimization process has stagnated
		exitFlag = 1;
		break
	elseif (iter > opt.maxIters)
		% max iterations reached
		exitFlag = 2;
		break
	end

	%% compute the gradient
	%===========================================================
	if opt.display; fprintf("calculating gradient: \n"); end
	%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	gradCost = cputime;
	gradTime = 0;
	tic;
	%...........................................................
	% backwards propagate
	phiBwd = propagate(phi,tgrid,field,propArgs{:},"direction","bwd","display",opt.display,"saveAll",true);
	% compute dP
	dP = dP_analytic(psiFwd,phiBwd,tgrid,mu,M,"propagator",opt.propagator);
	% apply pulse envelope
	dP_env = opt.envelope .* dP;
	%...........................................................
	gradTime = gradTime + toc;
	gradCost = cputime - gradCost;
	%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	% save
	if (saveProgress == true)
		saveStruct.dP(1,iter) = {dP};
		saveStruct.gradTime(1,iter) = {gradTime};
		saveStruct.gradCost(1,iter) = {gradCost};
	end
	
	% plot dP(t)
	if opt.plotting && isgraphics(ax3)
		plot(tNorm, dP_env, "Parent",ax3);
		drawnow;
	end

	% exit if max iterations reached
	if (iter > opt.maxIters)
		exitFlag = 2;
		break
	end

	%% compute the step-size "gamma"
	%===========================================================
	% estimate the ideal update step-size
	gamma = (sqrt(P) - P)/sum(abs(dP.*dP_env),"all");

	if opt.linesearch
		% employ a line-search to determine the optimal step-size
		%-----------------------------------------------------------
		gamma = opt.learnRate * gamma;
		if opt.display; fprintf("step-size (approx): %#11.4g \n",gamma); end
		%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		searchCost = cputime;
		searchTime = 0;
		tic;
		%...........................................................
		optFunc = @(x) calc_objective(field+(x.*dP_env));
		[gamma,~,gammaHist,valueHist] = line_search(optFunc,[0,gamma],[P,NaN], ...
			"bounds",[0,Inf],opt.searchOpts{:});
		%...........................................................
		searchTime = searchTime + toc;
		searchCost = cputime - searchCost;
		%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		if opt.display; fprintf("step-size (actual): %#11.4g \n",gamma); end
	else
		% use a static step-size
		gamma = opt.learnRate * gamma;
		if opt.display; fprintf("step-size: %#11.4g \n",gamma); end
	end
	
	% save
	if saveProgress
		saveStruct.gamma(1,iter) = {gamma};
		if opt.linesearch
			saveStruct.searchTime(1,iter) = {searchTime};
			saveStruct.searchCost(1,iter) = {searchCost};
			saveStruct.gammaHist(1,iter) = {gammaHist};
			saveStruct.valueHist(1,iter) = {valueHist};
		end
	end
	
	%% update the control field
	field = field + (gamma * dP_env);

	%% LOOP END
end

%% post process
%===========================================================
% save final information
if (saveProgress == true)
	saveStruct.iters = iter;
	saveStruct.exitFlag = exitFlag;
	saveStruct.optField = field;
	saveStruct.psiFwd = psiFwd;
end

% report exit condition
if opt.display
	switch exitFlag
		case 0; disp("optimization successful: overlap threshold reached.");
		case 1; disp("optimization stopped: optimizer reached local maximum.");
		case 2; disp("optimization stopped: maximum iterations reached.");
	end
end

%% internal functions
%***********************************************************
	function [P] = calc_objective(field)
	% forwards propagate
	psi_T = propagate(psi,tgrid,field,propArgs{:},"direction","fwd","display",false,"saveAll",false);
	
	% return the overlap
	P = overlap(phi,M,psi_T);
	end
%***********************************************************
%% DONE
end