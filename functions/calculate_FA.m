function [FA_values,DA_face,signed_diff,colormapData] = calculate_FA(all_faces,template,type)
%Inputes
% all_faces: contains the vertices of the already registered faces
% each row contains one face (vector of X, vector of Y and vector of Z values)
% template: Structure of vertices and faces of template mesh; m*3 array of xyz
% coordinates and their connectiones
% type: 'face'
%
%Outputs:
%FA_values: average fluctuating asymmetry per face
%DA_face: the directional asymmetry values per vertex
%signed_diff: signed difference between a face and its reflection,
%corrected for average DA
%colormapData: fluctuating asymmetry per vertex for each face ( each row
%contains the data for one face)


%   by Omid Ekrami, 2018

switch type
    case 'face'
        template.vertices = template.vertices- repmat(mean(template.vertices),size(template.vertices,1),1);
        tol = 0.0001;
        left=find(template.vertices(:,1)<-tol); %left side of the face
        midline = find(abs(template.vertices(:,1))<tol);
        pairs(1,:)=(1:length(left));
        for i=1:length(left)
            temp = template.vertices(i,:);
            temp(1) = -temp(1);
            ind(i) = find(ismembertol(template.vertices(length(left)+length(midline):end,:),temp,tol,'ByRows',1));
        end
        pairs(2,:)=ind+length(left)+length(midline)-1; % coantains the pair of corresponding vertices on the sides of the template face
    case 'body'
end

template.vertices = template.vertices - repmat(mean(template.vertices),size(template.vertices,1),1);

delta_ffactor = 1;
ffactor=0;
reference=template.vertices;
average=zeros(size(reference));

%% Weighted General Procrustes Analysis   
for num_face=1:size(all_faces,1)
        face = [all_faces(num_face,1:size(template.vertices,1));...
                all_faces(num_face,size(template.vertices,1)+1:2*size(template.vertices,1));...
                all_faces(num_face,2*size(template.vertices,1)+1:end)]';
        ref_face = get_reflection(face,pairs); %gives the mirror of the face with respect to x=0 plane
        average = average+face+ref_face;
        ref_faces(num_face,:) = reshape(ref_face,3*size(ref_face,1),1);
end
reference= average/(2*size(all_faces,1)); %average of all faces before Procrustes
pro_faces=all_faces;

while abs(delta_ffactor)>0.00001
    average=zeros(size(reference));
    for num_face=1:size(pro_faces,1)
        face = [pro_faces(num_face,1:size(template.vertices,1));...
                pro_faces(num_face,size(template.vertices,1)+1:2*size(template.vertices,1));...
                pro_faces(num_face,2*size(template.vertices,1)+1:end)]';
            %ref_face = get_reflection(face,pairs);  
        ref_face = [ref_faces(num_face,1:size(template.vertices,1));...
                    ref_faces(num_face,size(template.vertices,1)+1:2*size(template.vertices,1));...
                    ref_faces(num_face,2*size(template.vertices,1)+1:end)]';
        weights = ones(1,size(reference,1));
        change=1;
        while change>0.0001
            [~,face,~]=absor(face',reference','weights',weights); 
            %code can be found on Matlab file exchange under "Absolute Orientation - Horn's method by Matt J"
            [~,ref_face,~]=absor(ref_face',reference','weights',weights);
            face=face';
            ref_face=ref_face';
            Err= get_dist(face,ref_face);
            sigma=sqrt(sum(weights*Err)/sum(weights));
            lambda = (1/(2*pi*sigma))*exp(-2);
            new_weights = normpdf(Err,0,sigma)./(normpdf(Err,0,sigma)+lambda);
            change = mean(abs(weights-new_weights'));
            weights=new_weights';
        end
            average = average+ref_face+face;
            pro_faces(num_face,:) = reshape(face,3*size(face,1),1);
            ref_faces(num_face,:) = reshape(ref_face,3*size(face,1),1);
    end
    old_reference=reference;
    reference= average/(2*size(all_faces,1));
    old_ffactor = ffactor;
    ffactor=sum(get_dist(reference,old_reference));
    delta_ffactor = (ffactor-old_ffactor);
end

%% Calculate DA
switch type
    case 'face'
        colormapData = zeros(size(all_faces,1),size(all_faces,2)/3);
        signed_diff = zeros(size(all_faces));
        DA_face=0;
        for num_face=1:size(pro_faces,1)              
            pro_face = [pro_faces(num_face,1:size(template.vertices,1));...
                  pro_faces(num_face,size(template.vertices,1)+1:2*size(template.vertices,1));...
                  pro_faces(num_face,2*size(template.vertices,1)+1:end)]';
            ref_face = [ref_faces(num_face,1:size(template.vertices,1));...
                        ref_faces(num_face,size(template.vertices,1)+1:2*size(template.vertices,1));...
                        ref_faces(num_face,2*size(template.vertices,1)+1:end)]';
            DA_face =DA_face+(pro_face-ref_face);
        end
        DA_face =DA_face/size(pro_faces,1); % Directional Asymmetry
        
        for num_face=1:size(pro_faces,1)
            pro_face = [pro_faces(num_face,1:size(template.vertices,1));...
                  pro_faces(num_face,size(template.vertices,1)+1:2*size(template.vertices,1));...
                  pro_faces(num_face,2*size(template.vertices,1)+1:end)]';
            ref_face = [ref_faces(num_face,1:size(template.vertices,1));...
                        ref_faces(num_face,size(template.vertices,1)+1:2*size(template.vertices,1));...
                        ref_faces(num_face,2*size(template.vertices,1)+1:end)]';
            
            %difference between the face and its reflection        
            subtraction = sqrt(sum((pro_face-ref_face).^2,2));
           
            colormapData(num_face,:) = subtraction;
            FA_values(num_face) = mean(subtraction); %average FA per face
            signed_diff(num_face,:) =reshape(pro_face-ref_face,3*size(pro_face,1),1);

       end
case 'body'

end

FA_values=FA_values';
end
