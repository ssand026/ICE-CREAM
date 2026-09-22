function [TDM] = get_TDM(states, mu, opt)
% gives the transition dipole moment couplings for the provided field operators in mu
arguments
	states (:,:) double {mustBeFinite}
	mu     cell
	opt.noDiags (1,1) logical = false;
	opt.normalize (1,1) logical = false;
end
[numPts,numKet] = size(states);
numMu = numel(mu);

TDM = zeros(numKet,numKet,numMu);
for dd = 1:numMu
	% check the dipole-moment operator
	if any(size(mu{dd})~=numPts) && ~isscalar(mu{dd})
		error("ERROR: the dipole moment matrix mu{"+dd+"} has an invalid size")
	end
	
	tdm = (states' * mu{dd} * states);
	
	% normalize the tdm by the transition overlaps
	if opt.normalize
		row = sum(abs(tdm),1);
		col = sum(abs(tdm),2);
		tdm = 1/2 * (tdm./(~row + row) + tdm./(~col + col));
	end
	
	% ignore self-coupling
	if opt.noDiags
		tdm = tdm - diag(diag(tdm));
	end
	
	% append to output
	TDM(:,:,dd) = tdm;
end
% DONE
end