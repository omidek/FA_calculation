function mesh=clip_sphere(mesh, radius, centeroid)
%%cropping the mesh using a sphere
%   by Omid Ekrami, 2018

%find the points out of the sphere
    temp = repmat(centeroid,size(mesh.vertices,1) ,1);
    dist = get_dist(temp,mesh.vertices);
    outside = find(dist>radius);
    
%fix the indices in the faces
    reductions = zeros(size(mesh.faces));
    for i=1:length(outside)
        try
            temp_idx = find(mesh.faces>outside(i) & outside(i+1)>=mesh.faces); %fixes the other faces in the mesh
        catch
            temp_idx = find(mesh.faces>outside(i)); %fixes the other faces in the mesh
        end
        reductions(temp_idx) = reductions(temp_idx)-i;
    end
    extra_faces = ismember(mesh.faces,outside);
    extra_faces = sum(extra_faces,2);
    temp_idx = find(extra_faces~=0); %finds the mutual faces
    mesh.faces(temp_idx,:)=[]; %removes mutual faces
    reductions(temp_idx,:)=[];
    mesh.vertices(outside,:)=[]; %removes vertices outside the sphere
    mesh.faces = mesh.faces+reductions;
    unreferenced = find(~ismember(1:size(mesh.vertices,1),mesh.faces));
    reductions2 = zeros(size(mesh.faces));
    for i=1:length(unreferenced)
        try
            temp_idx = find(mesh.faces>unreferenced(i) & unreferenced(i+1)>=mesh.faces); %fixes the other faces in the mesh
        catch
            temp_idx = find(mesh.faces>unreferenced(i));
        end
        reductions2(temp_idx) = reductions2(temp_idx)-i;
    end
    mesh.vertices(unreferenced,:)=[]; %removes vertex2
    mesh.faces = mesh.faces+reductions2;
end