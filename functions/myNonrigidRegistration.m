function [registered]=myNonrigidRegistration(source,target,iterations,type)

% INPUT
% -target: Structure of vertices and faces of target mesh; n*3 array of xyz
% coordinates and their connectiones
% -source: Structure of vertices and faces of template mesh; m*3 array of xyz
% coordinates and their connectiones
% -iterations: number of iterations;
% -type: face, body or skull mesh

% OUTPUT
% -registered: registered template vertices on target mesh. Faces are not affected 

%   by Omid Ekrami, 2018

clf

%% finds the vertices on the edge of the mesh %% codes can be found on matlab file exchange
% Stable Sampling of Point Clouds for ICP Registration by Tolga Birdal
[source_edges]=detectedges(source.faces); 
[target_edges]=detectedges(target.faces);


%% plotting the meshes
trisurf(target.faces,target.vertices(:,1),target.vertices(:,2),target.vertices(:,3),'Facecolor','b','Edgecolor','none');
hold
light
set(gca, 'visible', 'off')
view(0,90)
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
alpha(0.4)
material dull
h=trisurf(source.faces,source.vertices(:,1),source.vertices(:,2),source.vertices(:,3),'FaceColor','y','Edgecolor','none');
pause (0.05)
%% General deformation

switch type
    case 'face'
    % parameters
    deg = (-6.5:4.4/iterations:-2.1); %temprature for deterministic annealing/ lower for all-to-one
    k = ceil(0.5.^(deg)); %number of corresponding vertices for each surface vertex (for KNN)
    
    maxNrPoints=500; 
    minNrPoints =25;
    nrPoints= ceil(maxNrPoints-1:(minNrPoints-maxNrPoints+1)/iterations:minNrPoints); %number of neighbouring points that affect each vertex

    outlier_dist=ones(1,iterations)*2; %distance used to select outliers
    outlier_angle = 0.9; % normal angle used for outlier detection
    outlier_weight=3;
    
    var1=30:-(30-1)/iterations:1;
    var2 = 2:(0.5-2)/iterations:0.5;
    
    case 'body'

end
%% finds surface normals
TR = triangulation(target.faces,target.vertices);
normalsT = vertexNormal(TR);

% Finding the maxNrPoints neighboing vertices to the vertex on source mesh
global_displacement = zeros(size(source.vertices));
original_source = source.vertices;
[indS_neighbors,Dsource]=knnsearch(source.vertices,source.vertices,'K',maxNrPoints);
dist_multiplier=1;

for counter=1:iterations
    
    TRS = triangulation(source.faces,source.vertices); 
    normalsS=vertexNormal(TRS); %surface normals for target surface
    
    %finding correspondences on target face for each vertex on source
    [indT_corr,Dtar_corr]=knnsearch([target.vertices,normalsT],[source.vertices,normalsS],'K',k(counter)); 
    
    init_targ=size(target.vertices,1);
    
    % Outlier detection, correction for the holes
    avg_dist = mean(mean(Dtar_corr));
    dist_SD = mean(std(Dtar_corr));
    orientation_weight = zeros(size(Dtar_corr));
    outlier=cell(k(counter),1);
    outlier2=cell(k(counter),1);
    outlier3=cell(k(counter),1);
    for j=1:k(counter)
        outlier{j}=find(ismember(indT_corr(:,j),target_edges)); %removes the correspondence with the edges of the mesh
        outlier2{j}=find(Dtar_corr(:,j)>(avg_dist+outlier_dist(counter)*dist_SD)); % removes correspondence with the vertices that are too far
        if exist('arm_ids')
            arm_outliers = find(ismember(outlier2{j},arm_ids));
            wrong_outliers = arm_outliers(find(Dtar_corr(outlier2{j}(arm_outliers),j)<(avg_dist+outlier_dist_arms(counter)*dist_SD)));
            outlier2{j}(wrong_outliers)=[];
        end
        orientation_weight(:,j) = (-acos(dot(normalsT(indT_corr(:,j),:),normalsS,2)))/pi+1;
        outlier3{j}=find(orientation_weight(:,j)<outlier_angle); %removes correspondence of faces with normal differences more than some angle
        orientation_weight(outlier3{j},j)=1;
        pp1=size(target.vertices,1);
        target.vertices=[target.vertices;source.vertices(outlier{j},:);source.vertices(outlier2{j},:);source.vertices(outlier3{j},:)];
        indT_corr(outlier{j},j)=pp1+(1:size(outlier{j},1))';
        indT_corr(outlier2{j},j)=pp1+(size(outlier{j},1))+(1:size(outlier2{j},1))';
        indT_corr(outlier3{j},j)=pp1+(size(outlier{j},1))+(size(outlier2{j},1))+(1:size(outlier3{j},1))';
        Dtar_corr(outlier{j},j)=0.01; %  the weight is set to zero
        Dtar_corr(outlier2{j},j)=outlier_weight;
        Dtar_corr(outlier3{j},j)=outlier_weight;
    end
    
    %finding the non overlapping areas of the surfaces
    non_overlap= find(sum((Dtar_corr-1000),2)==0);
    if ~isempty(non_overlap)
        Dtar_corr(non_overlap,:)=repmat(0.01,length(non_overlap),size(Dtar_corr,2));
    end
    
    % calculating the weigths of the K corresponding vertices on the target surface to the source surface based on distance
    correspondence_weights = (1/sqrt(2*pi*var2(counter)))*exp(-(1/(2*var2(counter)))*(dist_multiplier*Dtar_corr).^2);
    correspondence_weights = orientation_weight.*correspondence_weights;
    correspondence_weights = correspondence_weights./repmat(sum(correspondence_weights,2),1,size(correspondence_weights,2));
    Targettempset=0;
    
    for j=1:k(counter)
        Targettempset=Targettempset+horzcat(correspondence_weights(:,j).*target.vertices(indT_corr(:,j),1),...
                                            correspondence_weights(:,j).*target.vertices(indT_corr(:,j),2),...
                                            correspondence_weights(:,j).*target.vertices(indT_corr(:,j),3));
    end

    target.vertices=target.vertices(1:init_targ,:);

    % finding initial correspondence
    [indS_corr,Dsource_corr]=knnsearch([source.vertices,normalsS],[target.vertices,normalsT],'K',k(counter));
    init_source=size(source.vertices,1);
    
    % Outlier detection, correction for the holes
    avg_dist = mean(mean(Dsource_corr));
    dist_SD = mean(std(Dsource_corr));
    orientation_weight = zeros(size(Dsource_corr));
    
    for j=1:k(counter)
        outlier{j}=find(ismember(indS_corr(:,j),source_edges)); %removes the correspondence with the edges of the mesh
        outlier2{j}=find(Dsource_corr(:,j)>(avg_dist+outlier_dist(counter)*dist_SD)); % removes correspondence with the vertices that are too far
        orientation_weight(:,j) = (-acos(dot(normalsS(indS_corr(:,j),:),normalsT,2))/pi+1);
        outlier3{j}=find(orientation_weight(:,j)<outlier_angle); %removes correspondence of faces with normal differences more than some angle
        orientation_weight(outlier3{j},j)=1;
        pp1=size(source.vertices,1);
        source.vertices=[source.vertices;target.vertices(outlier{j},:);target.vertices(outlier2{j},:);target.vertices(outlier3{j},:)];
        indS_corr(outlier{j},j)=pp1+(1:size(outlier{j},1))';
        indS_corr(outlier2{j},j)=pp1+(size(outlier{j},1))+(1:size(outlier2{j},1))';
        indS_corr(outlier3{j},j)=pp1+(size(outlier{j},1))+(size(outlier2{j},1))+(1:size(outlier3{j},1))';
        Dsource_corr(outlier{j},j)=0.01; %  the weight is set to zero
        Dsource_corr(outlier2{j},j)=outlier_weight;
        Dsource_corr(outlier3{j},j)=outlier_weight;
    end
    
    %finding the non overlapping areas of the surfaces
    non_overlap= find(sum((Dsource_corr-1000),2)==0);
    if ~isempty(non_overlap)
        Dsource_corr(non_overlap,:)=repmat(0.01,length(non_overlap),size(Dsource_corr,2));
    end
    
    % calculating the weigths of the K corresponding vertices on the target surface to the source surface based on distance
    correspondence_weights = (1/sqrt(2*pi*var2(counter)))*exp(-(1/(2*var2(counter)))*(dist_multiplier*Dsource_corr).^2);
    correspondence_weights = orientation_weight.*correspondence_weights;
    correspondence_weights = correspondence_weights./repmat(sum(correspondence_weights,2),1,size(correspondence_weights,2));
    Sourcetempset=0;
    
    for j=1:k(counter)
        Sourcetempset=Sourcetempset+horzcat(correspondence_weights(:,j).*source.vertices(indS_corr(:,j),1),...
                                            correspondence_weights(:,j).*source.vertices(indS_corr(:,j),2),...
                                            correspondence_weights(:,j).*source.vertices(indS_corr(:,j),3));
    end

    source.vertices=source.vertices(1:init_source,:);

    
    % obtaining the displacement vectors
    starting = [source.vertices;Sourcetempset];
    finishing = [Targettempset;target.vertices];
    
    % calculating the weights of the neighboring vertices (the locality of the deformation)
    [indM_neighbors,Dmask]=knnsearch(starting,source.vertices,'K',nrPoints(counter));
    neighbor_weights = (1/sqrt(2*pi*var1(counter)))*exp(-(1/(2*var1(counter)))*(dist_multiplier*Dmask(:,1:nrPoints(counter))).^2);
    neighbor_weights = neighbor_weights./repmat(sum(neighbor_weights,2),1,size(neighbor_weights,2));
        
    fieldVector = finishing-starting;

    % Applying Gaussian on displacement vectors (viscousity)
    for i=1:size(source.vertices,1)
        for j=1:nrPoints(counter)
            deformationField(j,:)=neighbor_weights(i,j)*fieldVector(indM_neighbors(i,j),:);
        end
        global_displacement(i,:)=global_displacement(i,:)+sum(deformationField,1); %adding the effect of all the nearest points
    end
    %% Elasticity 
    neighbor_weights = (1/sqrt(2*pi*var1(counter)))*exp(-(1/(2*var1(counter)))*(dist_multiplier*Dsource(:,1:nrPoints(counter))).^2);
    neighbor_weights = neighbor_weights./repmat(sum(neighbor_weights,2),1,size(neighbor_weights,2));
    for i=1:size(source.vertices,1)
        for j=1:nrPoints(counter)
            sourceVtemp(j,:)=neighbor_weights(i,j)*global_displacement(indS_neighbors(i,j),:);
        end
        temp_GD(i,:)=sum(sourceVtemp,1); %adding the effect of all the nearest points
    end 
    global_displacement = temp_GD;
    source.vertices = original_source + global_displacement;
    if mean(get_dist(starting,finishing))<0.01
        break
    end
    %% displaying the mesh        
    delete(h)
    h=trisurf(source.faces,source.vertices(:,1),source.vertices(:,2),source.vertices(:,3),'FaceColor','y','edgealpha','0');   
    pause (0.05)
end
registered=source.vertices;
