function [faces] = boundary_faces(tri)
% Returns the set of faces that lie on the domain boundary
arguments
	tri (:,:) {mustBeInteger,mustBePositive}
end
[numVtx,~] = size(tri);

% get the set of all faces in the mesh
tri = sort(tri,1);
faceCombs = ncombsk(numVtx,numVtx-1);
faces = reshape(tri(faceCombs.',:),numVtx-1,[]).';

% return faces that only appear once in the triangulation as faces
% in the mesh interior will be found twice in the triangulation
isUnique = isunique(faces,1);
faces = faces(isUnique,:);
% DONE
end