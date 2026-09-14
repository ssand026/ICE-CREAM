function [newFaces,flipped] = align_normals(faces)
% ALIGN_NORMALS Aligns mesh normals to all point in a consistent direction.
%
% [newFaces,flipped] = ALIGN_NORMALS(faces) takes an [N x 3] array of surface
% triangles "faces" and returns a re-ordered set of faces "newFaces" such that
% all triangles have normals pointing in a similar direction. The second output
% "flipped" is an [N x 1] logical array denoting which faces were flipped in
% order to align the surface normals. This assumes face normals are calculated
% via the right-hand rule such that the vertex ordering matters.
arguments
	faces (:,:) double
end

% circumvent the costly overhead of "unique(...)" using:
%	[C,ia,ic] = matlab.internal.math.uniquehelper(A,doSort,isFirst,byRows)
fast_unique = @(x,opts) matlab.internal.math.uniquehelper(x,opts(1),opts(2),opts(3));

% check face list validity
[numFaces,numVtx] = size(faces);
if (numVtx ~= 3)
	error("ERROR: face list must have two or three vertices per face.")
end


% Get the list of all edges. Edges on the surface boundary will belong to
% only one face. To align the normals of the surface, each non-boundary edge
% should occur in the face list in ascending and descending order once, when 
% face normals are defined via the right-hand rule.
edges = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])];
faceIndex = repmat((1:numFaces).',3,1);
numEdges = size(edges,1);

% find edges where vertices are in ascending order
sortedEdges = (edges(:,1) < edges(:,2));

% get the list of unique edges, referenced by a group number
edges(~sortedEdges,:) = edges(~sortedEdges,[2 1]);
[~,~,edgeGroups] = fast_unique(edges,[false,true,true]);
numGroups = max(edgeGroups);

% find which edges are paired and already aligned
isPaired = (accumarray(edgeGroups,1,[numGroups,1])==2);
isAligned = (accumarray(edgeGroups,sortedEdges,[numGroups,1]) == 1);
isUnified = isPaired & isAligned;
isEdgeUnified = isUnified(edgeGroups);

% get the subset of paired edges and their location/face-index
pairSubset = isPaired(edgeGroups);
pairIndex = find(pairSubset);
pairedEdges = edgeGroups(pairSubset);
pairedFaces = faceIndex(pairSubset);

% for each paired edge, get the complementing face
[~,edgeOrder] = sort(pairedEdges);
setA = edgeOrder(1:2:end);
setB = edgeOrder(2:2:end);
edgeComplement = nan(numEdges,1);
edgeComplement(pairIndex(setA)) = pairedFaces(setB); 
edgeComplement(pairIndex(setB)) = pairedFaces(setA);

%% Collect connected faces/edges
% March from the first face to each of its nicely (asc/desc) connected
% neighbour faces. Label each connected "set" of faces.
faceSet = zeros(numFaces,1,"uint32");
lockedFaces = false(numFaces,1);
visitedEdges = false(numEdges,1);

currentSet = 0;
currFaces = [];
while any(~lockedFaces)
    % If we're not connected to anything, we must start a new set
    if isempty(currFaces)
        currFaces = find(~lockedFaces,1);
        currentSet = currentSet + 1;
    end
    lockedFaces(currFaces) = true;
    faceSet(currFaces) = currentSet;
    % Grab the edges of the current faces
    currEdgeIndx = bsxfun(@plus, currFaces, 0:numFaces:numEdges-1);
    % Find which edges are nicely connected to unvisited faces
    currEdgeIndxToFollowMask = isEdgeUnified(currEdgeIndx) & ~visitedEdges(currEdgeIndx);
    % Show that we've visited all edges of the current faces
    visitedEdges(currEdgeIndx) = true;
    % Determine the new faces we would reach if we stepped via nice edges
    linkedFaces = edgeComplement(currEdgeIndx(currEdgeIndxToFollowMask));
    currFaces = linkedFaces(~isnan(linkedFaces) & ~lockedFaces(linkedFaces));
end

%% Work out which sets need to be flipped and which stay the same
[setsTouched, setsToFlip] = deal(false(currentSet,1));
currentSets = [];
while any(~setsTouched)
    % If no current sets, pick the first one available. Any (next) sets
    % found touching it will need to be flipped
    if isempty(currentSets)
        currentSets = find(~setsTouched,1,'first');
        flipTheNextSet = true;
    end
    % We've now touched the current sets. Find these sets' faces.
    setsTouched(currentSets) = true;
    setsFaceInds = find(ismember(faceSet, currentSets));
    % Find edges that border these faces (includes edges from other sets)
    edgeIndsSharingBorder = find(ismember(edgeComplement, setsFaceInds));
    % Find all the faces that share these edges
    [faceNosSharingBorder,~] = ind2sub([numFaces,numVtx], edgeIndsSharingBorder);
    unqFaceNosSharingBorder = unique(faceNosSharingBorder);
    % Avoid flipping faces on any sets already touched
    unqFaceNosToFlipMask = ~setsTouched(faceSet(unqFaceNosSharingBorder));
    faceNosToFlip = unqFaceNosSharingBorder(unqFaceNosToFlipMask);
    % These will only be the border faces. Get the sets they belong to
    setNosToFlip = unique(faceSet(faceNosToFlip));
    % Flip those sets if we should. The first (root) set WON'T have been
    % flipped. Any touching it WILL get flipped. Any touching *those* will
    % already be in the same direction as the root set, so they WON'T be
    % flipped. The next WILL, WON'T, WILL, WON'T, etc.
    if flipTheNextSet
        setsToFlip(setNosToFlip) = true;
    end
    flipTheNextSet = ~flipTheNextSet;
    % Let the loop continue using the sets we just flipped as a source.
    currentSets = setNosToFlip;
end

%% Perform the actual flipping of all sets that require it
flipped = ismember(faceSet, find(setsToFlip));
newFaces = faces;
newFaces(flipped,:) = newFaces(flipped,[3 2 1]);

%% DONE
end