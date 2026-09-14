function [mu_new,field_scl] = rescale_dipoles(psi_a, psi_b, mu)
% Scales the dipole-moment operators so that the transition 
% has a dipole moment magnitude << 1
arguments
	psi_a (:,1) double
	psi_b (:,1) double
	mu    (:,:) cell
end
SCL = 0.1;

% find the dipole-moment magnitude
tdm = zeros(numel(mu),1);
for dd = 1:numel(mu)
	tdm(dd) = abs(psi_a' * mu{dd} * psi_b);
end
tdm = norm(tdm,2);

% In order for energy to be conserved:
%	(field_new * mu_new) = (field_old * mu_old)
% where:
%	mu_new = (SCL/tdm) * mu_old
% Therefore: 
%	field_new = (tdm/SCL) * field_old
field_scl = tdm/SCL;
mu_new = cell(size(mu));
for dd = 1:numel(mu)
	mu_new{dd} = (1/field_scl) * mu{dd};
end

% DONE
end