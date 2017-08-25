function [ H ] = est_homography(video_pts, logo_pts)
% est_homography estimates the homography to transform each of the
% video_pts into the logo_pts
% Inputs:
%     video_pts: a 4x2 matrix of corner points in the video
%     logo_pts: a 4x2 matrix of logo points that correspond to video_pts
% Outputs:
%     H: a 3x3 homography matrix such that logo_pts ~ H*video_pts
% Written for the University of Pennsylvania's Robotics:Perception course

% YOUR CODE HERE

A = zeros(8, 9);
for i = 1:4
    x1 = video_pts(i,1);
    x2 = video_pts(i,2);
    xp1 = logo_pts(i,1);
    xp2 = logo_pts(i,2);
    ax = [-x1, -x2, -1, 0, 0, 0, x1*xp1, x2*xp1, xp1];
    ay = [0, 0, 0, -x1, -x2, -1, x1*xp2, x2*xp2, xp2];

    A(i*2-1, :) = ax;
    A(i*2, :) = ay;
end

[U, S, V] = svd(A);

H = transpose(reshape(V(:,end), [3, 3]));

end

