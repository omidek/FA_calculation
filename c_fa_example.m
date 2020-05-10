%loading the mask and the set of sample faces
load('sym_KU.mat');
load('sample_faces.mat') % a set of 10 sample faces, synthetically made from a weighted average of real faces.
mask=sym;

%Calculating the signed difference between a face and its mirror
[~,~,signed_diff,~] = calculate_FA(sample_faces,mask,'face'); 

% the DA vector, calculated as the average of the signed differences
DA_vector = mean(signed_diff,1); 
DA_direction  = DA_vector'/norm(DA_vector);

% projection of the signed differences onto the DA direction and corrected for average DA, this produces F-DA score for each face
FDA_scores = signed_diff*DA_direction - norm(DA_vector);

% calculating the vector of projections for each face
DA_projections = repmat(FDA_scores,1,size(signed_diff,2)).*repmat(DA_direction',size(signed_diff,1),1);

%removing the DA projection for each face from its signed difference. This
%gives the asymmetry vector for each face, corrected for individual DA
%effects.
DA_residuals = signed_diff-DA_projections;

%calculates the average C-FA for each face
for i=1:size(DA_residuals,1)
    face = reshape(DA_residuals(i,:),size(DA_residuals,2)/3,3);
    C_FA(i) = mean(sqrt(sum(face.^2,2)));
end
