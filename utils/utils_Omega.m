function Omega = utils_Omega(w)
% UTILS_OMEGA Computes the 4x4 Omega matrix for quaternion kinematics
%
% Omega(w) is defined such that dq/dt = 0.5 * Omega(w) * q
%
% Reference: Lefferts, Markley, Shuster (1982), Eq. 16
%
% Inputs:
%   w - 3x1 angular velocity vector [wx; wy; wz]
%
% Outputs:
%   Omega - 4x4 matrix
%
% Notation:
%   q = [q1; q2; q3; q4] where q4 is the scalar part.
%   Omega = [ -[w x],  w
%             -w',     0 ]
%
% Coordinates:
%   matrix is:
%   [  0   w3 -w2  w1
%     -w3  0   w1  w2
%      w2 -w1  0   w3
%     -w1 -w2 -w3  0 ]

    wx = w(1);
    wy = w(2);
    wz = w(3);

    Omega = [  0,   wz, -wy,  wx;
              -wz,  0,   wx,  wy;
               wy, -wx,  0,   wz;
              -wx, -wy, -wz,   0 ];
end
