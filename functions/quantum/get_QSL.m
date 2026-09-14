function [QSL] = get_QSL(kind,states,H,M)
% returns the quantum speed limit (QSL) for transitions between every
% combination of states given in the basis "psi", given as an array of 
% column-vectors. The first input determines the type of QSL returned
% by the function, i.e. Mandelstam-Tamm, Margolus-Levitin, or Fermi.
arguments (Input)
	kind (1,1) {mustBeMember(kind,["tamm","levi","fermi"])}
	states (:,:) double
	H double {mustBeSquare}
	M double {mustBeSquare} = [];
end
arguments (Output)
	QSL (:,:) double
end

% normalize the state vectors
states = norm_wavefunc(states,M);

% compute the energy of each state
H_psi = (H * states);
expval_H = real(sum(conj(states) .* H_psi,1));

% compute the state overlap angles
if iseye(M)
	qAngle = abs(acos(abs(states' * (M * states))));
else
	qAngle = abs(acos(abs(states' * states)));
end

% get the quantum speed limit QSL(i,j) for transitions between psi(:,i)->psi(:,j)
switch kind
	case "tamm"
		% return the Mandelstam-Tamm QSL
		expval_HH = real(sum(conj(states) .* (H * H_psi),1));
		E_unc = sqrt(expval_HH - expval_H.^2); % energy uncertainty
		QSL = qAngle ./ abs(E_unc + E_unc.')/2;
	case "levi"
		% return the Margolus-Levitin QSL
		E_min = min(expval_H); % energy of the lowest-occupied state
		E_avg = (expval_H + expval_H.')/2; % average energy
		QSL = qAngle ./ abs(E_avg - E_min);
	case "fermi"
		% return the Heisenberg uncertainty QSL
		QSL = 4 * qAngle ./ abs(expval_H - expval_H.');
end

% fix invalid QSL values
QSL(~isfinite(QSL)) = NaN;

% DONE
end