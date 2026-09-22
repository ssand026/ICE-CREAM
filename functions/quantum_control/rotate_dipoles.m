function [mu_uv, T] = rotate_dipoles(psi_i, psi_f, mu_xy)
% Finds the linear combination of dipole-moments that has the strongest coupling
% to the given state-state transition.
arguments
	psi_i (:,1) double
	psi_f (:,1) double
	mu_xy (:,1) cell
end

% get direction of the strong-coupling axis
tdm = zeros(numel(mu_xy),1);
for ii = 1:numel(mu_xy)
	tdm(ii) = (psi_f' * (mu_xy{ii} * psi_i));
end
tdm = tdm / norm(tdm,2);

% apply phase shift to maximize real component
nRe = norm(real(tdm),2);
nIm = norm(imag(tdm),2);
phase = exp(-1i*angle(nRe + 1i*nIm));
tdm = phase .* tdm;

% extract strong axis from tdm
strongAxis = real(tdm);
strongAxis = strongAxis / norm(strongAxis,2);

% find orthogonal/weak axes
weakAxes = null(strongAxis.');

% generate the transformation matrix 'T' such that
%	[Fs; Fu; Fv] = T * [Fx; Fy; Fz]
T = [strongAxis, weakAxes].';
T = T ./ vecnorm(T,2,2);

% generate the new dipole moment operators
mu_uv = cell(size(mu_xy));
for ii = 1:numel(mu_xy)
	for jj = 1:numel(mu_xy)
		if isempty(mu_uv{ii})
			mu_uv{ii} = T(ii,jj) * mu_xy{jj};
		else
			mu_uv{ii} = mu_uv{ii} + T(ii,jj) * mu_xy{jj};
		end
	end
end
% DONE
end