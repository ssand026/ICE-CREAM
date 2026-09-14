function [xOpt, yOpt, xList, yList, exitFlag] = line_search(func, xPts, yPts, opt)
% LINE_SEARCH Golden-section search for a local minimum or maximum of the function.
%
% The algorithm combines golden-section interval refinement with
% automatic interval expansion when the current optimum lies at an
% interval boundary. Searches may therefore continue beyond the initial
% interval, subject to the specified search bounds.
%
% Syntax
% =======
% [xOpt,yOpt] = LINE_SEARCH(func,xPts) searches for a local maximum of
% the scalar function FUNC within the interval specified by xPts. xPts
% must contain two ascending values [xMin,xMax].
%
% [xOpt,yOpt] = LINE_SEARCH(func,xPts,yPts) supplies previously
% evaluated endpoint values yPts = [func(xMin), func(xMax)], avoiding
% redundant function evaluations.
%
% [xOpt,yOpt,xList,yList,exitFlag] = LINE_SEARCH(...) additionally
% returns the complete search history and termination status.
%
% Input Arguments
% ================
% func - Function handle returning a real-valued scalar.
% xPts - Initial search interval [xMin,xMax].
% yPts - Optional endpoint values corresponding to xPts. Use NaN to force evaluation.
%
% Name-Value Arguments
% =====================
%  bounds - Hard search limits. The search interval may expand until these bounds are reached.
%  goal   - Optimization objective:
%              "max"  Search for a maximum (default)
%              "min"  Search for a minimum
%  display   - Set to true to display iteration progress.
%  parallel  - Evaluate independent function values using PARFOR when a parallel pool is available.
%  threshold - Relative improvement threshold used to declare convergence.
%  maxStall  - Maximum number of consecutive iterations without improvement before terminating.
%  minIters  - Minimum number of iterations before convergence tests are permitted.
%  maxIters  - Maximum number of search iterations.
%
% Output Arguments
% =================
% xOpt     - Location of the optimal point found.
% yOpt     - Function value at xOpt.
% xList    - List of points sampled during the linesearch
% yList    - List of the corresponding function values.
% exitFlag - Termination status:
%               0:  Convergence threshold reached
%               1:  Maximum iterations reached
%               2:  Maximum stall count reached
%
arguments
	func function_handle
	xPts (1,2) double {mustBeReal,mustBeFinite}
	yPts (1,2) double {mustBeReal} = [NaN,NaN];
	% optional arguments
	opt.bounds (1,2) double {mustBeReal} = [NaN,NaN];
	opt.goal (1,1) {mustBeMember(opt.goal,["min","max"])} = "max";
	opt.display (1,1) logical = true;
	opt.parallel (1,1) logical = true;
	% optimization parameters
	opt.threshold (1,1) double {mustBeInRange(opt.threshold,0,1)} = 0.05;
	opt.maxStall (1,1) uint16 {mustBePositive} = 5;
	opt.minIters (1,1) uint16 {mustBePositive} = 3;
	opt.maxIters (1,1) uint16 {mustBePositive} = 100;
end

%% initialize
%===========================================================
% ensure that the search-interval contains ascending values
if (xPts(1) >= xPts(2))
	error("ERROR: the initial search-interval must contain ascending values")
end

% ensure the search bounds capture the initial interval
xBnd = [min(opt.bounds(1),xPts(1)), max(opt.bounds(2),xPts(2))];

% set search objective
switch opt.goal
	case "min"; best = @(varargin) min(varargin{:}); dir = +1;
	case "max"; best = @(varargin) max(varargin{:}); dir = -1;
end

% initialize display format
if opt.display
	num2str = @(val) join(arrayfun(@(x)sprintf("%#10.4g",x),val)," ");
	fprintf("\n");
	fprintf("      x1        x2         x3         x4     |");
	fprintf("    f(x1)     f(x2)      f(x3)      f(x4)   \n");
	fprintf(" ============================================|");
	fprintf("=========================================== \n");
end

% check parallelism
if (opt.parallel==true) && isempty(gcp("nocreate"))
	opt.parallel = false;
end

% ensure minimimum iterations are completed
opt.maxIters = max(opt.minIters,opt.maxIters);

%% evaluate the initial search-interval points/values
%===========================================================
igr = (sqrt(5) - 1)/2; % the inverted golden-ratio: ~0.6180

% set initial x-values
x = [xPts(1), NaN, NaN, xPts(2)];
x(2) =  (igr)  * x(1) + (1-igr) * x(4);
x(3) = (1-igr) * x(1) +  (igr)  * x(4);

% define the validated function
f = @(x) fcheck(func(x));

% check which x-values require evaluation
y = [yPts(1), NaN, NaN, yPts(2)];
toEval = find(isnan(y));

if opt.parallel
	% evaluate y-values in parallel
	parfor ii = toEval
		y(ii) = f(x(ii));
	end
else
	% evaluate y-values sequentially
	for ii = toEval
		y(ii) = f(x(ii));
	end
end

% display progress
if opt.display
	fprintf(" %s ~%s \n",num2str(x),num2str(y))
end

% save data to output list
xList = x(:);
yList = y(:);

%% update search-interval until a local min/max is found
%===========================================================
stallSteps = 0;
prevBest = best(y);
convgLvl = NaN;
exitFlag = 1;

for iter = 1:opt.maxIters
	% choose a new point to evaluate the function at depending
	% on which sub-interval contains a local minima/maxima
	%-----------------------------------------------------------
	[~,best_x] = best(y);
	if (best_x == 1) && (x(1) > xBnd(1))
		% add new point before x(1)
		symbol = "<";
		xNew = x(1) - (x(2) - x(1))/igr;
		xNew = max(xNew,xBnd(1));
		yNew = f(xNew);
		x = [xNew, x(1), x(2), x(3)];
		y = [yNew, y(1), y(2), y(3)];

	elseif (best_x == 4) && (x(4) < xBnd(2))
		% add new point after x(4)
		symbol = ">";
		xNew = x(4) + (x(4) - x(3))/igr;
		xNew = min(xNew,xBnd(2));
		yNew = f(xNew);
		x = [x(2), x(3), x(4), xNew];
		y = [y(2), y(3), y(4), yNew];
		
	elseif (best_x == 1) || (best_x == 2)
		% add new point between x(1) and x(3)
		symbol = "-";
		xNew = (igr) * x(1) + (1-igr) * x(3);
		yNew = f(xNew);
		x = [x(1), xNew, x(2), x(3)];
		y = [y(1), yNew, y(2), y(3)];
		
	elseif (best_x == 3) || (best_x == 4)
		% add new point between x(2) and x(4)
		symbol = "+";
		xNew = (1-igr) * x(2) + (igr) * x(4);
		yNew = f(xNew);
		x = [x(2), x(3), xNew, x(4)];
		y = [y(2), y(3), yNew, y(4)];
	end
	
	% ensure x-values are in ascending order (disordering can
	% occur if the interval expansion hits a limit)
	[~,order] = sort(x);
	x = x(order);
	y = y(order);
	
	% save the new results
	xList(end+1,:) = xNew; %#ok
	yList(end+1,:) = yNew; %#ok
	
	% display progress
	if opt.display
		fprintf(" %s "+symbol+"%s \n",num2str(x),num2str(y));
	end


	% update convergence tests if the y-value has improved
	%-----------------------------------------------------------
	currBest = best(y);
	if (currBest ~= prevBest) && (best([currBest,prevBest]) == currBest)
		% calculate the relative improvement
		convgLvl = abs(currBest - prevBest)/(abs(currBest)+abs(prevBest));
		% update the previous best
		prevBest = currBest;
		stallSteps = 0;
	else
		stallSteps = stallSteps + 1;
	end
	
	% check all exit conditions
	%-----------------------------------------------------------
	isExpanding = ((symbol=="<") | (symbol==">")); % true if search interval was expanded
	isConverged = (convgLvl < opt.threshold); % true if the optimal y-value has converged
	isStagnant = (stallSteps >= opt.maxStall); % true if the optimal y-value stopped improving
	canExit = (~isExpanding && (iter > opt.minIters)); % true if minimum exit requirements are met
	
	if canExit && isConverged
		% exit if converged
		exitFlag = 0;
		break
	elseif canExit && (iter < opt.maxIters) && isStagnant
		% exit if convergence has stagnated
		exitFlag = 2;
		break
	end
end

%% post-process
%===========================================================
% interpolate test-points to find precise local minimum/maximum
searchOpt.Display = 'off';
searchOpt.TolX = 0;
searchOpt.TolFunc = 1e-12;
xOpt = fminbnd(@(xq) dir * interp1(xList,yList,xq,'spline'),min(xList),max(xList),searchOpt);
yOpt = f(xOpt);

% display the exit reason and optimized results
if opt.display
	switch exitFlag
		case 0; fprintf("\n convergence threshold reached: \n");
		case 1; fprintf("\n maximum iterations reached: \n");
		case 2; fprintf("\n maximum improvement attempts reached: \n");
	end
	fprintf("     x  = %s \n", strtrim(num2str(xOpt)));
	fprintf("   f(x) = %s \n", strtrim(num2str(yOpt)));
end

%% DONE
end

% external functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fx] = fcheck(fx)
% checks if the function evaluation fx = func(x) is valid
if ~isnumeric(fx) || ~isscalar(fx) || isnan(fx) || ~isreal(fx)
	error("ERROR: the function must return a real-valued numeric scalar")
else
	fx = double(fx);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%