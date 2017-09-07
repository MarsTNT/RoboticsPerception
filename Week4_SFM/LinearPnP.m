function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

N = size(x,1);

X = [X ones(N,1)];
x = [x ones(N,1)];

A = [];
for i=1:N
    skew = Vec2Skew(inv(K)*x(i,:)');
    a = [X(i,:)         zeros(1,4)     zeros(1,4);
         zeros(1,4)     X(i,:)         zeros(1,4);
         zeros(1,4)     zeros(1,4)     X(i,:)]; 
    A = [A; skew*a];
end

[U, S, V] = svd(A);

P = reshape(V(:,end), [4, 3])';
R_prime = P(:,1:3);
t = P(:,4);

[U, S, V] = svd(R_prime);

if det(U*V') > 0
    Rc = U*V';
    tc = t/S(1,1);
else
    Rc = -U*V';
    tc = -t/S(1,1);   
end

C = -Rc'*tc;
R = Rc;

end






