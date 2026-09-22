function [faces,regions] = region_bounds(tri,reg)
% Finds all surfaces that form the boundary between different regions of the mesh
%
% SEE ALSO: BOUNDARY_FACES
arguments
	tri (:,:) {mustBeInteger,mustBePositive}
	reg (1,:) {mustBeInteger,mustBeNonnegative}
end

% get the set of exterior faces
bndFaces = boundary_faces(tri);

if isscalar(reg) || all(reg==reg(1))
	% --- the mesh only contains a single region ---
	faces{1} = bndFaces;
	regions{1} = [0, reg(1)];
else
	% --- find faces that are shared amongst two regions  ---
	regVals = unique(reg(:));
	numReg = numel(regVals);

	% find the boundary faces for each region in the mesh
	facesInReg = cell(1,numReg);
	for ii = 1:numReg
		reg_ii = (reg == regVals(ii));
		facesInReg{ii} = boundary_faces(tri(:,reg_ii));
	end

	% find the faces shared between multiple regions
	sharedFaces = cell(numReg,numReg);
	sharedRegions = cell(numReg,numReg);
	for ii = 1:numReg
		for jj = ii:numReg
			if ii==jj
				% surface is on the mesh exterior
				sharedFaces{ii,jj} = intersect(facesInReg{ii},bndFaces,"rows");
				sharedRegions{ii,jj} = [0,regVals(ii)];
			else
				% surface is on the mesh interior
				sharedFaces{ii,jj} = intersect(facesInReg{ii},facesInReg{jj},"rows");
				sharedRegions{ii,jj} = [regVals(ii),regVals(jj)];
			end
		end
	end
	
	% return any surfaces with non-zero overlap
	hasOverlap = ~cellfun(@isempty,sharedFaces);
	faces = sharedFaces(hasOverlap);
	regions = sharedRegions(hasOverlap);
end
% DONE
end