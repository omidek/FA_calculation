function [registered,varargout] = initial_aligning (mask,target,flag, type, varargin)
% Aligning the two surfaces using the provided positioning landmarks
%   by Omid Ekrami, 2018

switch type
    case 'face'
    mask_landmarks = [1388 5779 3591 2248,4919];
    case 'body'

end

if nargin ==5
    target_landmarks = varargin{1};
elseif nargin == 6
    mask_landmarks = varargin{1};
    target_landmarks = varargin{2};
elseif nargin == 7
    mask_landmarks = varargin{1};
    target_landmarks = varargin{2};
    check_landmarks = varargin{3};
end


if flag==1

%% "Least-Squares Fitting of Two 3-D Point Sets" by K.S. Arun, T.S. Huang, and S.D. Blostein
    q1 = mask.vertices(mask_landmarks,:)';
    q2 = target.vertices(target_landmarks,:);

    R = absoluteOrientation(q1,q2);
    T = mean(target.vertices)' - R*mean(mask.vertices)';
    T_mat = [R,T;0 0 0 1];
    varargout{2} = T_mat;
    mask.vertices = T_mat*[mask.vertices';ones(1,size(mask.vertices,1))];
    mask.vertices = mask.vertices(1:3,:)';
    diff = target.vertices(target_landmarks(3),:)-mask.vertices(mask_landmarks(3),:);
    if exist('face_flag','var')
        mask.vertices = mask.vertices+repmat(diff,size(mask.vertices,1),1);
    end
    
    if exist('check_landmarks','var')
        num_lms = size(check_landmarks,1);
        temp_landmarks = [check_landmarks(:,1:3);check_landmarks(:,4:6);check_landmarks(:,7:9);check_landmarks(:,10:12)];
        temp_landmarks = T_mat*[temp_landmarks';ones(1,size(temp_landmarks,1))];
        temp_landmarks = temp_landmarks(1:3,:)';
        check_landmarks = [temp_landmarks(1:num_lms,:),temp_landmarks(num_lms+1:2*num_lms,:),temp_landmarks(2*num_lms+1:3*num_lms,:),temp_landmarks(3*num_lms+1:4*num_lms,:)];
        if exist('face_flag','var')
            check_landmarks = check_landmarks+repmat(diff,size(check_landmarks,1),4);
        end
    end
end


%% procrustes registeration

TR = triangulation(target.faces,target.vertices); 
normalsT = vertexNormal(TR);
TRS = triangulation(mask.faces,mask.vertices); 
normalsS=vertexNormal(TRS);

[IDX1(:,1),IDX1(:,2)]=knnsearch([target.vertices,normalsT],[mask.vertices,normalsS]);
[IDX2(:,1),IDX2(:,2)]=knnsearch([mask.vertices,normalsS],[target.vertices,normalsT]);
IDX1(:,3)=1:length(mask.vertices(:,1));
IDX2(:,3)=1:length(target.vertices(:,1));

m1=mean(IDX1(:,2));
s1=std(IDX1(:,2));
IDX1=IDX1(IDX1(:,2)<(m1+1.67*s1),:); % 97.5 percentile point of the normal distribution

m1=mean(IDX2(:,2));
s1=std(IDX2(:,2));
IDX2=IDX2(IDX2(:,2)<(m1+1.67*s1),:); % 97.5 percentile point of the normal distribution

Datasetsource=vertcat(mask.vertices(IDX1(:,3),:),mask.vertices(IDX2(:,1),:));
Datasettarget=vertcat(target.vertices(IDX1(:,1),:),target.vertices(IDX2(:,3),:));

[~,~,transform] = procrustes(Datasettarget, Datasetsource, 'reflection', 0,'scaling',0);
mask.vertices=transform.b*mask.vertices*transform.T+repmat(transform.c(1,:),size(mask.vertices,1),1);
varargout{3} = transform; 
if exist('check_landmarks','var')
    temp_landmarks = [check_landmarks(:,1:3);check_landmarks(:,4:6);check_landmarks(:,7:9);check_landmarks(:,10:12)];
    temp_landmarks = transform.b*temp_landmarks*transform.T+transform.c(1:size(temp_landmarks,1),:);
    varargout{1} = [temp_landmarks(1:num_lms,:),temp_landmarks(num_lms+1:2*num_lms,:),temp_landmarks(2*num_lms+1:3*num_lms,:),temp_landmarks(3*num_lms+1:4*num_lms,:)];
end

registered = mask.vertices;
