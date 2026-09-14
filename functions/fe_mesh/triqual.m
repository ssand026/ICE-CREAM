function [aspect] = triqual(pts,tri)
% Returns the edge-ratio based quality measure for each simplex in the mesh
combo = ncombsk(size(tri,1),2);
edg_len = zeros(size(combo,1),size(tri,2));
for ii = 1:size(combo,1)
	node_a = pts(:,tri(combo(ii,1),:));
	node_b = pts(:,tri(combo(ii,2),:));
	edg_len(ii,:) = sqrt(sum((node_a - node_b).^2,1));
end
aspect = (min(edg_len,[],1) ./ max(edg_len,[],1));
end