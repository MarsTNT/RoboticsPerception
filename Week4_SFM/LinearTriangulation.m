function X = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)
%% LinearTriangulation
% Find 3D positions of the point correspondences using the relative
% position of one camera from another
% Inputs:
%     C1 - size (3 x 1) translation of the first camera pose
%     R1 - size (3 x 3) rotation of the first camera pose
%     C2 - size (3 x 1) translation of the second camera
%     R2 - size (3 x 3) rotation of the second camera pose
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Outputs: 
%     X - size (N x 3) matrix whos rows represent the 3D triangulated
%       points

N = size(x1,1);

T1 = -R1 * C1;
T2 = -R2 * C2;
P1 = K * [R1 T1];
P2 = K * [R2 T2];

X = zeros(N,3);

for i = 1:N  
    x1_h = [x1(i,:) 1];
    x2_h = [x2(i,:) 1];
    A = [Vec2Skew(x1_h)*P1;
         Vec2Skew(x2_h)*P2];
    
    [U,S,V] = svd(A);
    X(i, :) = V(1:3,end)/V(end,end);
end


end