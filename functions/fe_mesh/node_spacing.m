function [d] = node_spacing(n,vol,area,len,efficiency)
% Estimates the average node-separation "d" for a region containing "n" nodes. 
% Nodes can lie on the region boundary or in the region's interior. Using the
% following syntax will return the node separation for:
% 
%	NODE_SPACING(n,a,b,[])
%		3D region with volume "a" and surface-area "b"
%	
%	NODE_SPACING(n,a,[],[])
%		spherical region with volume "a" 
%
%	NODE_SPACING(n,[],b,c)
%		2D region with area "b" and perimiter length "c"
%
%	NODE_SPACING(n,[],b,[])
%		circular region with area "b"
%	
%	NODE_SPACING(n,[],[],c)
%		1D region with length "c"
%
% The optional argument "efficiency" sets how close the packing is to optimal,
% where a value of 1 assumes optimal circle/sphere packing. Defaults to 0.9 for
% 2D systems and 0.75 for 3D systems.
arguments
	n (1,1) double {mustBeGreaterThanOrEqual(n,2)}
	vol  double {mustBeScalarOrEmpty,mustBeReal,mustBeNonnegative,mustBeFinite} = [];
	area double {mustBeScalarOrEmpty,mustBeReal,mustBeNonnegative,mustBeFinite} = [];
	len  double {mustBeScalarOrEmpty,mustBeReal,mustBeNonnegative,mustBeFinite} = [];
	efficiency double {mustBeScalarOrEmpty,mustBeInRange(efficiency,0,1)} = [];
end

if ~isempty(vol) && isempty(len)
	% region is 3D
	%-----------------------------------------------------------
	if isempty(efficiency); efficiency = 0.75; end
	closePack = pi/(3*sqrt(2)); % optimal sphere-packing ratio in 3D
	fillRatio = efficiency * closePack;

	% surface area if region is spherical
	areaMin = 4*pi* (vol/(4/3*pi))^(2/3);
	area = max([area,areaMin]);
	
	% As nodes can lie on the region boundary, the effective sphere-packing
	% volume is approximately equal to:
	%	effectiveVol = (vol + area * d/2)
	% If we define the relation between node-separation and filled volume:
	%	(fillRatio * effectiveVol)/n = 4*pi/3*(d/2)^3
	% We obtain:
	%	(d^3) - (6/pi*fillRatio/n * area/2)*d - (6/pi*fillRatio/n * vol) = 0
	beta = (6/pi) * fillRatio/n;
	d = roots([1, 0, -beta*area/2, -beta*vol]);
	d = max(real(d));

elseif ~isempty(area) && isempty(vol)
	% region is 2D
	%-----------------------------------------------------------
	if isempty(efficiency); efficiency = 0.9; end
	closePack = pi/(2*sqrt(3)); % optimal circle-packing ratio in 2D
	fillRatio = efficiency * closePack;

	% perimeter length if region is circular
	lenMin = 2*pi*sqrt(area/pi);
	len = max([len,lenMin]);

	% As nodes can lie on the region boundary, the effective circle-packing
	% area is approximately equal to:
	%	effectiveArea = (area + len * d/2)
	% If we define the relation between node-separation and filled area:
	%	(fillRatio * effectiveArea)/n = pi*(d/2)^2 
	% We obtain:
	%	(d^2) - (4/pi*fillRatio/n * len/2)*d - (4/pi*fillRatio/n * area) = 0
	beta = (4/pi) * fillRatio/n;
	d = roots([1, -beta*len/2, -beta*area]);
	d = max(d);

elseif ~isempty(len) && isempty(vol) && isempty(area)
	% region is 1D
	%-----------------------------------------------------------
	d = len/(n-1);
else
	error("ERROR: invalid specification for the volume, area, and perimeter.")
end
% DONE
end