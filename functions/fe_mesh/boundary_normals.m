function [nvec,faces] = boundary_normals(pts,tri,faces)
% Computes the normal vector for each face on the domain boundary
% The resulting normals point outwards from the domain by default.
% Specifying a third argument only computes the normals for faces in the list.
%
% SEE ALSO: BOUNDARY_FACES, FACE_NORMALS
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	faces (:,:) {mustBeInteger,mustBePositive} = boundary_faces(tri);
end

[numDim,~] = size(pts);
[numVtx,numTri] = size(tri);
[numFace,~] = size(faces);

% check the validity of the face list
if nargin >= 3
	if size(faces,2)~=numDim
		error("ERROR: the number of vertices per face must equal the mesh dimensionality")
	end
	faces = sort(faces,2);
	if any(~isunique(faces,1))
		error("ERROR: the face list contains duplicates")
	end
end

% extract the portion of triangulation containing at least one face
tri = sort(tri,1);
combs = ncombsk(numVtx,numDim);
hasFaces = false(numTri,1);
for ii = 1:size(combs,1)
	subset_ii = tri(combs(ii,:),:).';
	hasFaces = hasFaces | ismember(subset_ii,faces,"rows");
end
tri = tri(:,hasFaces);

% check if boundary exists
if isempty(tri) || isempty(faces)
	% exit early
	nvec = [];
	faces = [];
	return
end

% find the location of each face in the triangulation and get the normal vectors
nvec = zeros(numDim,numFace);
inTri = false(numFace,1);
for ii = 1:size(combs,1)
	% check if each face is contained in the current permutation of vertices   
	set_ii = permute(tri(combs(ii,:),:),[2,1]);
	[isFace,faceIndex] = ismember(set_ii,faces,"rows");
	faceIndex = faceIndex(isFace,:);
		
	% get face positions
	inSet_ii = false(numFace,1);
	inSet_ii(faceIndex) = true;

	if ~any(isFace)
		% proceed if the current vertex permutation does not contain boundary faces
		continue
	elseif (nargin >= 3) 
		% check if the triangulation contains duplicate faces
		if any(~isunique(faceIndex)) || any(inTri & inSet_ii)
			error("ERROR: triangulation contains duplicates of one or more faces")
		end
	end
	
	% update match list
	inTri = (inTri | inSet_ii);

	% get the valid faces and their opposing nodes
	face_ii = set_ii(isFace,:);
	
	% get the nodes opposite of each face
	not_ii = setdiff([1:numVtx],combs(ii,:));
	not_ii = tri(not_ii,isFace).';
	
	% compute normals and orient outwards from the domain
	norm_ii = face_normals(pts,face_ii);
	side_ii = pts(:,face_ii(:,1))-pts(:,not_ii);
	norm_ii = sign(dot(norm_ii,side_ii,1)).*norm_ii;
	
	nvec(:,inSet_ii) = norm_ii;
end

% check if any faces were omitted
if any(~inTri)
	warning("some of the listed faces were missing from the triangulation")
end
% DONE
end