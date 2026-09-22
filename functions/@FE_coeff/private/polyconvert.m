function [vtxToPoly,polyToVtx] = polyconvert(vtxPoly)
% returns the matrices that convert between the physical values of a function at
% the sub-simplex midpoints and the weights of a set of barycentric polynomials
% where each row of "vtxPoly" denotes the powers of the barycentric coordinates
% that define one basis function in the basis
%
% Suppose we have a simplex defined by the set of vertices {xi, xj, ...}.
% If the function "f[x]" is evaluated at the simplex's vertices, as well as the 
% midpoints of each sub-simplex such that:
%    f_i = f[xi],  f_ij = f[(xi + xj)/2],   f_ijk = f[(xi + xj + xk)/3], 
% we can define f[x] using the barycentric coordinates {λi[x], λj[x], ...}:
%    f[x] = [λi[x] ... λi[x]*λj[x] ...] * [b_i ... b_ij ... ].'
%
% This function evaluates the matrix "vtxToPoly" and "polyToVtx" such that:
%  [b_i ... b_ij ... b_ijk ...].' = vtxToPoly*[f_i ... f_ij ... f_ijk ...].'
%  [f_i ... f_ij ... f_ijk ...].' = polyToVtx*[b_i ... b_ij ... b_ijk ...].'
%
% as coefficients are typically evaluated/plotted using the physical values "f"
% wheras while the barycentric weights "b" are required by many other 
% calculations, such as integration/differentiation
arguments
	vtxPoly (:,:) double
end

% get the matrix P, where P_jk = 1 if the basis function in the j-th row of 
% "vtxPoly" depends on the k-th basis function. P will be lower-triangular if 
% the rows of "vtxPoly" contain basis-functions of increasing order
numPoly = size(vtxPoly,1);
P = zeros(numPoly,numPoly);
for ii = 1:numPoly
	P(ii,:) = ~any(~vtxPoly(ii,:) & vtxPoly,2);
	P(ii,ii) = false;
end

% the function value 'f' at the nth-degree simplex-point will equal the sum of
% all the basis-function weights 'b', times the value of the basis function at
% the simplex-point: (1/n)^(deg). Thus to get the conversion from 'f' to 'b' we
% rearrange: 
%   f_i = (b_i)
%   f_ij = (b_ij)*(1/2)^2 + (b_i+b_j)*(1/2)
%   f_ijk = (b_ijk)*(1/3)^3 + (b_ij+b_jk+b_ik)*(1/3)^2 + (b_i+b_j+b_k)*(1/3)
% so that
%   b_i = f_i
%   b_ij = (2^2)*(f_ij) - (2^1)*(b_i+b_j)
%   b_ijk = (3^3)*(f_ijk) - (3^1)*(b_ij+b_jk+b_ik) - (3^2)*(b_i+b_j+b_k)
% etc. where the previous values of b_{n} are used to determine b_{n+1}
deg = sum(vtxPoly,2);
vtxToPoly = eye(numPoly,numPoly);
for nn = 1:max(deg)
	vtxToPoly(deg==nn,:) = -((nn.^(nn-deg.')) .* P(deg==nn,:)) * vtxToPoly;
	vtxToPoly = vtxToPoly + (nn^nn)*diag(deg==nn);
end

% invert to get the back-conversion matrix
polyToVtx = inv(vtxToPoly);
% DONE
end