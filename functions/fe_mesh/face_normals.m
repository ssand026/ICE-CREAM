function [nvec] = face_normals(pts,faces)
% Computes the non-oriented normal vector for each face in the list
% The output is a [numDim x numFace] array
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	faces (:,:) {mustBeInteger,mustBePositive}
end
[numDim,~] = size(pts);
[numFace,~] = size(faces);
if size(faces,2)~=numDim
	error("ERROR: the number of face vertices must equal the mesh dimensionality")
end

% get the side-length vectors
svec = zeros(numDim,numDim-1,numFace);
for ii = 1:numFace
	svec(:,:,ii) = pts(:,faces(ii,1:end-1)) - pts(:,faces(ii,end));
end

% normalize and append with ones
nS = vecnorm(svec,2,1);
svec = svec ./ (nS + ~nS);
svec(:,end+1,:) = sqrt(numDim);

% get normal vector (via Gram-Schmidt orthogonalization)
Q = zeros(size(svec));
for jj = 1:size(Q,2)
	Q(:,jj,:) = svec(:,jj,:);
	if jj > 1
		proj = sum(Q(:,1:jj-1,:).*Q(:,jj,:),1);
		proj = sum(Q(:,1:jj-1,:).*proj,2);
		Q(:,jj,:) = Q(:,jj,:) - proj;
	end
	% normalize
	nQ = vecnorm(Q(:,jj,:),2,1);
	Q(:,jj,:) = Q(:,jj,:) ./ (nQ + ~nQ);
end
nvec = permute(Q(:,end,:),[1,3,2]); % numDim x numFace array
% DONE
end