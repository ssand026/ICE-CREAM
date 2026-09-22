function [pulse] = pulse_envelope(tgrid,type,params)
% Returns a pulse envelope of the specified type using the given parameters
% The returned pulses have a maximum value of one
%
% pulse types:
%	gauss: gaussian pulse using the given pulse midpoint and standard deviation
%	sech: hyperbolic secant pulse using the given pulse midpoint and fwhm
%	sine: pulse with a sin(pi/2*t)^2 distribution over the given ramp-up/down times
%	trap: trapezoidal pulse with a linear distribution during the ramp-up/down times
arguments
	tgrid (1,:) double {mustBeReal,mustBeFinite}
	type (1,1) {mustBeMember(type,["sine","gauss","trap","sech"])}
	params (1,2) double {mustBeReal,mustBeFinite}
end

% check that tgrid contains ascending values
if any(diff(tgrid,1,2) <= 0)
	error("ERROR: time-grid must be an array of ascending values.")
end

% remove initial shift
tgrid = tgrid - tgrid(1);
T = tgrid(end);

% calculate the ramping function start/end points
if matches(type,["sine","trap"])
	tL = params(1); % ramp-up time
	tR = params(2); % ramp-down time
	if (tL + tR) > T
		% rescale ramp times to fit the total duration
		scl = T/(tL+tR);
		tL = scl * tL;
		tR = scl * tR;
	end
	
	% ramp-up endpoint
	indxL = find(tgrid < tL,1,"last");
	if isempty(indxL); indxL = 0; end
	rampL = tgrid(1:indxL)/tL;
	
	% ramp-down startpoint
	indxR = find(tgrid > (T-tR),1,"first");
	if isempty(indxR); indxR = numel(tgrid)+1; end
	rampR = (T-tgrid(indxR:end))/tR;
end

% return the pulse envelope
switch type
	case "gauss"
		tMid = params(1); % pulse midpoint
		tDev = params(2); % pulse standard deviation
		pulse = exp(-(tgrid-tMid).^2 /(2*tDev.^2));
	case "sech"
		tMid = params(1); % pulse midpoint
		tWid = params(2); % pulse FWHM
		pulse = sech((tgrid-tMid) * asech(1/2)/tWid);
	case "sine"
		pulse = [sinpi(rampL/2).^2, ones(1,indxR-indxL-1), sinpi(rampR/2).^2];
	case "trap"
		pulse = [rampL, ones(1,indxR-indxL-1), rampR];
end
% DONE
end