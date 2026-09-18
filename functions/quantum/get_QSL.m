function [QSL] = get_QSL(kind,psi,H,M)
% GET_QSL returns the quantum speed limits for transitions between the basis states "psi".
% 
% Returns the an array containing the quantum speed-limits (QSLs) for transitions between each
% combination of states in the basis set "psi", which is given as an array of column vectors. The
% first input determines the whether the returned QSL array is obtained from the Mandelstam-Tamm or 
% from the Margolus-Levitin QSL. Specify the overlap matrix M if the basis is non-orthogonal.
arguments (Input)
	kind (1,1) {mustBeMember(kind,["tamm","levi"])}
	psi (:,:) double
	H double {mustBeSquare}
	M double {mustBeSquare} = [];
end
arguments (Output)
	QSL (:,:) double {mustBeReal}
end

% normalize the state vectors
psi = norm_wavefunc(psi,M);

% compute the state overlap angles
if iseye(M)
	qAngle = acos(abs(psi' * psi));
else
	qAngle = acos(abs(psi' * (M * psi)));
end
qAngle = abs(qAngle);

% compute the energy of each state
H_psi = (H * psi);
expval_H = sum(conj(psi) .* H_psi,1);

% get the quantum speed limit QSL(i,j) for transitions between psi(:,i)->psi(:,j)
switch kind
	case "tamm"
		% return the Mandelstam-Tamm QSL
		if iseye(M)
			expval_HH = real(sum(conj(psi) .* (H * H_psi),1));
		else
			expval_HH = real(sum(conj(psi) .* (H * (M\H_psi)),1));
		end
		E_unc = sqrt(abs(expval_HH - expval_H.^2)); % energy uncertainty
		QSL = 2*qAngle ./ (E_unc + E_unc.');
	case "levi"
		% return the Margolus-Levitin QSL
		QSL = 2*qAngle ./ abs(expval_H - expval_H.');
end

% remove invalid QSL values
QSL(~isfinite(QSL)) = NaN;

% DONE
end
