function Xi = utils_Xi(q)
% UTILS_XI Computes the 4x3 Xi matrix
%
% This matrix relates angular velocity to quaternion rates or errors.
% Relationships:
%   dq/dt = 0.5 * Xi(q) * w
%   Omega(w) * q = Xi(q) * w
%
% Inputs:
%   q - 4x1 quaternion vector [q1; q2; q3; q4] (scalar last)
%
% Outputs:
%   Xi - 4x3 matrix
%
% Formula:
%   Xi(q) = [ q4*eye(3) + [q1:3 x]
%             -q1:3' ]

    q1 = q(1); 
    q2 = q(2); 
    q3 = q(3); 
    q4 = q(4);

    % Cross product matrix of vector part [q1:3 x]
    % [ 0  -q3  q2
    %   q3  0  -q1
    %  -q2  q1  0 ]
    qv_cross = [ 0,  -q3,  q2;
                 q3,  0,  -q1;
                -q2,  q1,  0 ];

    Xi = [ q4 * eye(3) + qv_cross;
           -q(1:3)' ];
end
