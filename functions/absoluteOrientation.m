function R = absoluteOrientation(q1,q2)

%finds the rotation matrix to transfer q1 on to q2
%q1 3Xn
%q2 nX3

meanLeft = mean(q1,2);
meanRight = mean(q2,1)';

 M = q1*q2;
 
 M = M - size(q1,2)*meanLeft*meanRight';

delta = [M(2,3) - M(3,2); M(3,1) - M(1,3);M(1,2) - M(2,1)];
 
 N = [trace(M) delta'; delta (M+M'-trace(M)*eye(3))];
 
[eigenVectors, eigenValues] = eig(N);  
[~,index] = max(diag(eigenValues));
  
R = quat2rotm(eigenVectors(:,index)');

translation = meanRight - R*meanLeft;
T = [R, translation; [0, 0, 0, 1]];

end
