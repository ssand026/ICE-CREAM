function [statFunc] = stats_func(data,dims)
% STATS_FUNC generates a function that evaluates various statistical measures for the given data.
% 
% Calculates the minimum, maximum, magnitude, median, mean, and standard-deviation of the
% data along the specified dims, and returns a function that evaluates any string-expression
% containing the labels {MIN,MAX,MAG,MED,AVG,STD}. If dims is not specified, these quantities are
% evaluated along the first non-singulardimension of data.
arguments
	data double 
	dims (1,:) double {mustBeInteger,mustBePositive} = [];
end

if (numel(data) <= 1)
	error("ERROR: data must contain multiple elements.")
end

if isempty(dims)
	dims = find(size(data)~=1,1,"first");
end

% get statistics for nodeData
dMin = min(data,[],dims);
dMax = max(data,[],dims);
dMag = max(abs(dMin),abs(dMax));
dAvg = mean(data,dims);
dMed = median(data,dims);
dStd = std(data,0,dims) + eps(dMag);

% return evaluator function
statFunc = @(str) feval(str2func("@(MIN,MAX,MAG,MED,AVG,STD)"+str),dMin,dMax,dMag,dMed,dAvg,dStd);
end