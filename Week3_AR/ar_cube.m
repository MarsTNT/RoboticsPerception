function [proj_points, t, R] = ar_cube(H,render_points,K)
%% ar_cube
% Estimate your position and orientation with respect to a set of 4 points on the ground
% Inputs:
%    H - the computed homography from the corners in the image
%    render_points - size (N x 3) matrix of world points to project
%    K - size (3 x 3) calibration matrix for the camera
% Outputs: 
%    proj_points - size (N x 2) matrix of the projected points in pixel
%      coordinates
%    t - size (3 x 1) vector of the translation of the transformation
%    R - size (3 x 3) matrix of the rotation of the transformation
% Written by Stephen Phillips for the Coursera Robotics:Perception course

% YOUR CODE HERE: Extract the pose from the homography
[t, R] = compute_pose(H);
if t(3) < 0
	disp(['t3 (z) is negative, so computing alternative pose']);
	[t, R] = compute_pose(-1*H);
end

% YOUR CODE HERE: Project the points using the pose
[N,M] = size(render_points);
proj_points = zeros(N,2);

for j = 1:N
	X = [render_points(j,1); render_points(j,2); render_points(j,3)];
	Xc = K*(R*X + t);
	proj_points(j,1) = Xc(1) / Xc(3);
	proj_points(j,2) = Xc(2) / Xc(3);
end
   
end


function [t, R] = compute_pose(H)
h1 = H(:,1);
h2 = H(:,2);
h3 = cross(h1,h2);
Rp = [ h1 h2 h3];

[U,S,V] = svd(Rp);
S = [ 1 0 0 ; 0 1 0; 0 0 det(U*V')];
R = U*S*V';

t = H(:,3)/norm(h1);
end
   