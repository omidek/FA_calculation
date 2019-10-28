function mesh=clip_sphere(mesh, radius, centeroid)
%%cropping the mesh using a sphere

mesh_out = mesh; % copy all other fields of mesh to the output mesh_out
v = mesh.vertices;
f = mesh.faces;

%% find the index of rows of the vertices to be deleted
temp = repmat(centeroid,size(v,1) ,1);
dist = get_dist(temp,v);
outside = find(dist>radius);
logicalRemoveVertices = ismember(1:length(v),outside); % logical vector: true for rows to be removed

%%  Create new vertice reference tags
vnew = v;
tagsOld = 1:length(v);
tagsNew = tagsOld; 
newCount = cumsum(~logicalRemoveVertices); % counts the vertices that are to remain
tagsNew(~logicalRemoveVertices) = newCount(~logicalRemoveVertices); % the newCount are the new reference tags
tagsNew(logicalRemoveVertices) = nan; 

%% Delete vertices
vnew(logicalRemoveVertices,:) = []; 

%% Delete faces
fnew = f; 
[IndexFacesRowDelete,~] = find(ismember(f,outside)); % finds the row index of the faces to delete
fnew(IndexFacesRowDelete,:) = []; % deletes faces referencing the removed vertices

%% Renumber faces
fnewRenumbered = tagsNew(fnew); 

%% New Patch
mesh_out.faces = fnewRenumbered;
mesh_out.vertices = vnew;

%% Recursively call to Remove additional Vertices that are unreferenced by the removal of faces
unreferenced = 1:length(vnew);
unusedVerticeList = unreferenced(~ismember(unreferenced,fnewRenumbered));
if ~isempty(unusedVerticeList)
    mesh_out = removeVerticesPatch(mesh_out,mesh_out.vertices(unusedVerticeList,:));
end
