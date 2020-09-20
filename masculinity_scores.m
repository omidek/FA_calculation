
% registered_scans is a n-by-m matrix, containing the standardized data of 
% n faces, each represented by m quasi-landmarks. By standardized
% we mean a facial mask was non-rigidly registered onto each of these
% faces. "Ekrami et al.,2018"
% Each row contains landmark data for one face
% [x1,x2,x3,...,xn,y1,y2,y3...,yn,z1,z2,z3...,zn]

%gender is a vector of n-by-1 containing the the sex of each face (eg.
%male=1 female=2)

% we use a combination of the first 10 pls component as they explain more
% than 85% of the covariance between the matrix of faces
% and the vector of genders
nComp = 10;
[PLSeigenvectors,eigenvectors2,scores,YS,BETA,PCTVAR,MSE,stats] = plsregress(registered_scans,gender,nComp);

%% calculating scores
% we only used the first pls component to calculate the masculinity scores
% as it contains 
comp_num=(1:10); 
masculinity_vector = PLSeigenvectors(:,comp_num)/norm(PLSeigenvectors(:,comp_num));
masc_scores = sum((registered_scans*masculinity_vector),2);
masc_scores=(masc_scores)/max(abs(masc_scores));

% male_ind and female_ind contain the indices of male and female faces
masc_scores = masc_scores-(mean(masc_scores(female_ind))+mean(masc_scores(male_ind)))/2;
