function [vtxData] = poly_to_vertex(baryData,poly)
% Given a field defined via the coefficients of a barycentric polynomial,
% returns the field values at the corresponding polynomial points (i.e. the
% midpoints of the simplex faces).
% Required as the polynomial coefficients are required for most computations
% involving FE scalar/vector/tensor fields, whereas the actual field values
% are necessary for evaluation and plotting.

% apply the transformation to the data
[~,B2F] = polyconvert(poly);
if isequal(B2F,eye(size(B2F)))
	% conversion matrix is the identity matrix
	vtxData = baryData;
else
	% apply conversion matrix
	vtxData = tensorprod(B2F,baryData,[2],[3]);
	vtxData = permute(vtxData,[2,3,1,4]);
end
% DONE
end