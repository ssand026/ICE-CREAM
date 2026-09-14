function [baryData] = vertex_to_poly(vtxData,poly)
% Given a field defined via its value on the vertices (and face-midpoints)
% of each simplex, return the coefficients of a barycentric polynomial, that
% fits the field values.
% Required as the polynomial coefficients are required for most computations
% involving FE scalar/vector/tensor fields, whereas the actual field values
% are necessary for evaluation and plotting.

% apply the transformation to the data
[F2B,~] = polyconvert(poly);
if isequal(F2B,eye(size(F2B)))
	% conversion matrix is the identity matrix
	baryData = vtxData;
else
	% apply conversion matrix
	baryData = tensorprod(F2B,vtxData,[2],[3]);
	baryData = permute(baryData,[2,3,1,4]);
end
% DONE
end