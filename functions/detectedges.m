function[edge_idx]=detectedges(faces)

v1 = faces(:,1);
v2 = faces(:,2);
v3 = faces(:,3);

edge1=sort([v1 v2]')';
edge2=sort([v1 v3]')';
edge3=sort([v2 v3]')';

%non_shared edges
edges=[edge1 ;edge2 ;edge3];
[~,ia]=unique(edges,'rows','stable');
non_shared=edges(ia,:);

%shared edges
shared=removerows(edges,ia);

borders = setdiff(non_shared,shared,'rows');
edge_idx=unique(reshape(borders,size(borders,1)*2,1));

