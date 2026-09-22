function [out] = data_contours(data,numShells,spacing,opt)
% DATA_CONTOURS Determine contours according to the slope of the sorted data distribution.
% 
% In regions where the data distribution is relatively flat, contour surfaces do not convey much
% information. In regions where the data distribution changes rapidly, contour surfaces will enclose
% a fairly small proportion of the data.
arguments
	data (1,:) double
	numShells (1,1) double {mustBeInteger,mustBePositive} 
	spacing (1,1) {mustBeMember(spacing,["volume","values"])}
	opt.minSlope (1,1) double {mustBeInRange(opt.minSlope,0,1)} = 1/2;
	opt.maxSlope (1,1) double {mustBeGreaterThan(opt.maxSlope,1)} = 10;
end

len = numel(data);
sdata = sort(data);

% choose bin edgess so that each bin contains a similar number of data points
nBins = ceil(sqrt(len));
y = unique(interp1(1:len,sdata,linspace(1,len,nBins)));
x = histcounts(sdata,"BinEdges",y,"Normalization","probability");
x = cumsum([0,x]);

% get the normalized slope of the data distribution
xMid = (x(1:end-1)+x(2:end))/2;
dydx = diff(y)./diff(x);
dydx = dydx/(sdata(end)-sdata(1));

% get left/right bounds of the valid distribution regions
valid = isbetween(dydx,opt.minSlope,opt.maxSlope);
indxL = find(diff([0,valid])==+1);
indxR = find(diff([valid,0])==-1);
[indxL,indxR] = deal(indxL(indxL~=indxR),indxR(indxL~=indxR));

% compute the contour placements
%-----------------------------------------------------------
switch spacing
	case "volume"
		% contours enclose equal proportions of the data
		xL = xMid(indxL);
		xR = xMid(indxR);
		numDivs = ((xR-xL)/sum(xR-xL) * (numShells+1));
		numDivs = round(numDivs);
		
		surfX = [];
		for ii = 1:numel(numDivs)
			xGrid = ((1:numDivs(ii))-1/2)/numDivs(ii);
			surfX = [surfX, xL(ii) + (xR(ii) - xL(ii))*xGrid]; %#ok
		end
		out = interp1(linspace(0,1,len),sdata,surfX);
	case "values"
		% contours are equally spaced among the data values
		yL = y(indxL);
		yR = y(indxR);
		numDivs = ((yR-yL)/sum(yR-yL) * (numShells+1));
		numDivs = round(numDivs);
		
		surfY = [];
		for ii = 1:numel(numDivs)
			yGrid = ((1:numDivs(ii))-1/2)/numDivs(ii);
			surfY = [surfY, yL(ii) + (yR(ii) - yL(ii))*yGrid]; %#ok
		end
		out = surfY;
end
% DONE
end