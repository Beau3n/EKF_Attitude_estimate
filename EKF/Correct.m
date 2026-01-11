function [x_corr, P_corr] = Correct(x_pred, P_pred, meas_vecs, ref_vecs, R)
% CORRECT Performs the measurement update (correction) step of the EKF
%
% Inputs:
%   x_pred     - Predicted state [q; b] (7x1)
%   P_pred     - Predicted error covariance (6x6)
%   meas_vecs  - Measured vectors in body frame (3xM)
%   ref_vecs   - Reference vectors in inertial frame (3xM)
%   R          - Measurement noise covariance for a single vector (3x3)
%                (Assumed same for all M measurements)
%
% Outputs:
%   x_corr     - Corrected state [q; b] (7x1)
%   P_corr     - Corrected error covariance (6x6)

    % Number of measurements
    M = size(meas_vecs, 2);
    
    if M == 0
        x_corr = x_pred;
        P_corr = P_pred;
        return;
    end

    q_pred = x_pred(1:4);
    b_pred = x_pred(5:7);
    
    % Compute Attitude Matrix from predicted quaternion
    A_pred = utils_AttitudeMatrix(q_pred);
    
    % Initialize aggregated Jacobian and Residual
    % For M measurements, z is 3*M x 1
    H = zeros(3*M, 6);
    res = zeros(3*M, 1);
    
    % Construct Big R matrix
    % Method 1: Block diagonal if R is 3x3
    % R_all = kron(eye(M), R); 
    % However, if M is large, this is inefficient. We process sequentially or block-wise.
    % Only 2-3 stars typically, so block diagonal is fine.
    
    for i = 1:M
        v_meas = meas_vecs(:, i);
        v_ref  = ref_vecs(:, i);
        
        % Predict measurement: z_hat = A * r
        v_hat = A_pred * v_ref;
        
        % Residual: z - z_hat
        % Note: Careful with normalization? 
        % Star tracker vectors are usually unit vectors. 
        % We assume inputs are normalized or raw measurements match the model.
        r_i = v_meas - v_hat;
        
        % Jacobian block for this measurement
        % H_i = [ [v_hat x], 0_3x3 ]
        % Defined such that z - z_hat = H * delta_x + v
        % z approx z_hat + [z_hat x] * delta_theta
        
        % Cross product matrix of v_hat
        vx = [ 0, -v_hat(3), v_hat(2);
               v_hat(3), 0, -v_hat(1);
              -v_hat(2), v_hat(1), 0 ];
              
        H_i = [vx, zeros(3, 3)];
        
        % Fill into large arrays
        row_start = (i-1)*3 + 1;
        row_end   = i*3;
        
        H(row_start:row_end, :) = H_i;
        res(row_start:row_end) = r_i;
    end
    
    % Measurement Noise Covariance
    % Assuming R is 3x3 for each vector, independent
    R_all = kron(eye(M), R);
    
    % Kalman Gain
    % K = P * H' * inv(H * P * H' + R)
    S = H * P_pred * H' + R_all;
    K = P_pred * H' / S; % / S is clearer than inv(S)
    
    % State Update (Error State)
    dx = K * res;
    
    dtheta = dx(1:3);
    dbias  = dx(4:6);
    
    % Update State
    % 1. Bias update: b = b + db
    b_corr = b_pred + dbias;
    
    % 2. Quaternion update
    % q (+) = q (-) + 0.5 * Xi(q) * dtheta
    % Then normalize
    Xi = utils_Xi(q_pred);
    q_corr = q_pred + 0.5 * Xi * dtheta;
    q_corr = q_corr / norm(q_corr);
    
    % Reassemble
    x_corr = [q_corr; b_corr];
    
    % Covariance Update
    % P = (I - K*H) * P
    % Joseph form is more numerically stable: P = (I-KH)P(I-KH)' + KRK'
    I = eye(6);
    P_corr = (I - K*H) * P_pred * (I - K*H)' + K * R_all * K';
    
end
