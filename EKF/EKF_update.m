function [x_hist, P_hist] = EKF_update(t, gyro_meas, star_meas, star_ref, x_init, P_init, Q, R)
% EKF_UPDATE 运行姿态 EKF 主循环
% Runs the main Attitude EKF loop
%
% Note: using the name 'filter' causes shadowing
%
% Inputs:
%   t          - 时间向量 Time vector (1xN)
%   gyro_meas  - 陀螺仪测量 Gyro measurements (3xN) [rad/s]
%   star_meas  - 星敏感器测量 Star measurements (3xMxN) [unit vectors in body]
%   star_ref   - 参考星向量 Reference star vectors (3xM) [unit vectors in inertial]
%   x_init     - 初始状态 Initial state [q; b] (7x1)
%   P_init     - 初始协方差 Initial covariance (6x6)
%   Q          - 过程噪声协方差 Process noise covariance (6x6)
%   R          - 测量噪声协方差 Measurement noise covariance (3x3 for single vector)
%
% Outputs:
%   x_hist     - 状态历史 State history (7xN)
%   P_hist     - 协方差历史 Covariance history (6x6xN)

    N = length(t);
    M = size(star_ref, 2); % Number of stars
    
    % Allocate memory
    x_hist = zeros(7, N);
    P_hist = zeros(6, 6, N);
    
    % Initialize
    x_curr = x_init;
    P_curr = P_init;
    
    % Time step (assume constant for simplicity, or calculate per step)
    dt_const = t(2) - t(1);
    
    % Progress bar
    % h_wait = waitbar(0, 'Running EKF...');
    
    for k = 1:N
        % 1. Time Update (Prediction)
        % Get current gyro measurement
        w_meas = gyro_meas(:, k);
        
        % Calculate dt (if variable)
        if k < N
            dt = t(k+1) - t(k);
        else
            dt = dt_const;
        end
        
        % Call Predict
        [x_pred, P_pred] = EKF_prediction(x_curr, P_curr, w_meas, dt, Q);
        
        % 2. Measurement Update (Correction)
        % Get current star measurements
        % star_meas is 3 x M x N
        z_k = star_meas(:, :, k);
        
        % Call Correct
        % Note: ref_vecs are constant in inertial frame
        [x_curr, P_curr] = Correct(x_pred, P_pred, z_k, star_ref, R);
        
        % 3. Store results
        x_hist(:, k) = x_curr;
        P_hist(:, :, k) = P_curr;
        
        % if mod(k, 100) == 0
        %    waitbar(k/N, h_wait);
        % end
    end
    
    % close(h_wait);

end
