function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x - size (N × 2) matrices whose row represents correspondence
%                      between the first, second, and third images where 
%                      N is the number of correspondences.
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 

N = size(X0,1);
X = zeros(N,3);
for i = 1:N
    X(i,:) = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), X0(i,:));
end

end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)

i = 0;
while i < 1
    J1 = Jacobian_Triangulation(C1, R1, K, X0);
    J2 = Jacobian_Triangulation(C2, R2, K, X0);
    J3 = Jacobian_Triangulation(C3, R3, K, X0);
    J = [J1' J2' J3']';

    b = [x1(1) x1(2) x2(1) x2(2) x3(1) x3(2)]';

    f1 = K*R1*(X0' - C1);
    f2 = K*R2*(X0' - C2);
    f3 = K*R3*(X0' - C3);
    
    fx = [f1(1)/f1(3) f1(2)/f1(3) f2(1)/f2(3) f2(2)/f2(3) f3(1)/f3(3) f3(2)/f3(3)]';
    
    dX = J'*J \ J'*(b - fx);
    
    X0 = X0 + dX';
    
    i = i +1;
end
X = X0;

end

function J = Jacobian_Triangulation(C, R, K, X)
f = K(1,1);
px = K(1,3);
py = K(2,3);
dudX = [f*R(1,1) + px*R(3,1), f*R(1,2) + px*R(3,2), f*R(1,3) + px*R(3,3) ];
dvdX = [f*R(2,1) + py*R(3,1), f*R(2,2) + py*R(3,2), f*R(2,3) + py*R(3,3) ];
dwdX = [R(3,1), R(3,2), R(3,3) ];
x = K*R*(X' - C);
J = [(x(3).*dudX - x(1).*dwdX)./(x(3)*x(3));
     (x(3).*dvdX - x(2).*dwdX)./(x(3)*x(3))];
end
