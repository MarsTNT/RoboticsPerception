function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

N = size(x1,1);

x1n = x1';
x2n = x2';

% Build the constraint matrix
A = [x2n(1,:)'.*x1n(1,:)'   x2n(1,:)'.*x1n(2,:)'  x2n(1,:)' ...
     x2n(2,:)'.*x1n(1,:)'   x2n(2,:)'.*x1n(2,:)'  x2n(2,:)' ...
     x1n(1,:)'              x1n(2,:)'             ones(N,1) ];  
 
[U, S, V] = svd(A);

% Fundamental matrix is the column of V corresponding to
% smallest singular value.
F = transpose(reshape(V(:,end), [3, 3]));

% Apply rank 2 constraint on F
[U, S, V] = svd(F);
F = U*diag([S(1,1) S(2,2) 0])*V';

end