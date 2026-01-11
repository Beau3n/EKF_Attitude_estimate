function A = utils_AttitudeMatrix(q)
% UTILS_ATTITUDEMATRIX Computes the Direction Cosine Matrix A from Quaternion q
%
% Transforms vectors from Inertial to Body frame: v_b = A * v_i
%
% Inputs:
%   q - 4x1 quaternion vector [q1; q2; q3; q4] (scalar last)
%
% Outputs:
%   A - 3x3 Attitude Matrix (DCM)
%
% Formula (Standard for q representing I->B rotation):
%   A = (q4^2 - |v|^2)I + 2*v*v' - 2*q4*[v x]

    v = q(1:3);
    s = q(4);
    
    v_norm_sq = v' * v;
    s_sq = s * s;
    
    % Cross product matrix [v x]
    vx = [ 0,    -v(3),  v(2);
           v(3),  0,    -v(1);
          -v(2),  v(1),  0 ];
          
    A = (s_sq - v_norm_sq) * eye(3) + 2 * (v * v') - 2 * s * vx;
    
end
