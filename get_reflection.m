function ref_face = get_reflection(vertices,pairs)
%Inputs:
%vertices: vertices of the original face
%pairs: a 2*n matrix containing the indices of the paired points

%outputs
%ref_face: the vertices of the reflected and relabled face

%   by Omid Ekrami, 2018
%% This piece of code finds the pairs in a symmetrical template (before mapping)
% template.vertices = template.vertices- repmat(mean(template.vertices),size(template.vertices,1),1);
% tol = 0.0001;
% left=find(template.vertices(:,1)<-tol); %left side of the face
% midline = find(abs(template.vertices(:,1))<tol);
% pairs(1,:)=(1:length(left));
% for i=1:length(left)
%     temp = template.vertices(i,:);
%     temp(1) = -temp(1);
%     ind(i) = find(ismembertol(template.vertices(length(left)+length(midline):end,:),temp,tol,'ByRows',1));
% end
% pairs(2,:)=ind+length(left)+length(midline)-1;
%%
    ref_face= vertices;
    ref_face(:,1)=-ref_face(:,1); %inverting the x values
    
    %relabling (swapping the points with their pairs
    temp = ref_face;
    ref_face(pairs(1,:),:)= temp(pairs(2,:),:);
    ref_face(pairs(2,:),:)= temp(pairs(1,:),:);
end
