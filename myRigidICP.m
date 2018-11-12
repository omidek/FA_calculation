function [registered] = myRigidICP(source,target,num_samples,scaleFlag)
%   Rigid Iterative Closest Points algorithm
%   by Omid Ekrami, 2018

    %% finds the vertices on the edge of the mesh
    [source_edges]=detectedges(source.faces); 
    [target_edges]=detectedges(target.faces);
    
    TR = triangulation(target.faces,target.vertices); 
    normalsT = vertexNormal(TR);
    TRS = triangulation(source.faces,source.vertices); 
    normalsS=vertexNormal(TRS);
    
    %% sample points selection %% codes can be found on matlab file exchange
    % Stable Sampling of Point Clouds for ICP Registration by Tolga Birdal
    [~,~,indS]=sample_pc_stable(source.vertices, normalsS, num_samples);
    [~,~,indT]=sample_pc_stable(target.vertices, normalsT, num_samples);
%%
    targetSample=target.vertices(indT,:);
    targetSampleNormals=normalsT(indT,:);

    change=20;
    threshold=0.001;
    while change>threshold
        sourceSample = source.vertices(indS,:);
        sourceSampleNormals=normalsS(indS,:);

        IDX1=[];
        IDX2=[];
        [IDX1(:,1),IDX1(:,2)]=knnsearch([target.vertices,normalsT],[sourceSample,sourceSampleNormals]);
        [IDX2(:,1),IDX2(:,2)]=knnsearch([source.vertices,normalsS],[targetSample,targetSampleNormals]);
        IDX1(:,3)=1:length(sourceSample);
        IDX2(:,3)=1:length(targetSample);

        [row,~] = find(ismember(IDX1(:,1),target_edges));
        IDX1(row,:)=[];
        [row,~]=find(ismember(IDX2(:,1),source_edges));
        IDX2(row,:)=[];

        Avg=mean(IDX1(:,2));
        SD=std(IDX1(:,2));
        IDX1=IDX1(IDX1(:,2)<(Avg+1.67*SD),:); % 67 percentile point of the normal distribution
        
        Avg=mean(IDX2(:,2));
        SD=std(IDX2(:,2));
        IDX2=IDX2(IDX2(:,2)<(Avg+1.67*SD),:); % 67 percentile point of the normal distribution

        Datasetsource=vertcat(sourceSample(IDX1(:,3),:),source.vertices(IDX2(:,1),:));
        Datasettarget=vertcat(target.vertices(IDX1(:,1),:),targetSample(IDX2(:,3),:));

        [~,~,transform] = procrustes(Datasettarget, Datasetsource, 'reflection', 0,'scaling',scaleFlag);
        old_source = source.vertices;
        source.vertices=transform.b*source.vertices*transform.T+repmat(transform.c(1,:),size(source.vertices,1),1);
        change =  mean(get_dist(source.vertices,old_source));

    end
    registered= source.vertices;
end