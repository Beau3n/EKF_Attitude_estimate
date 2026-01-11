function [x_pred, P_pred] = Predict(x_est, P_est, omega_meas, dt, Q)
% PREDICT 执行 EKF 的时间更新 (预测) 步骤
% Performs the time update (propagation) step of the EKF
%
% 算法基于 Lefferts, Markley, Shuster (1982)
% Algorithm based on Lefferts, Markley, Shuster (1982)
%
% Inputs:
%   x_est      - 当前状态估计 [q; b] (7x1) Current state estimate
%                q: 四元数 quaternion [q1, q2, q3, q4]'
%                b: 陀螺零偏 gyro bias [bx, by, bz]'
%   P_est      - 当前误差协方差 (6x6) Current error covariance
%   omega_meas - 陀螺仪测量的角速度 (3x1) Measured angular velocity from gyro
%   dt         - 时间步长 (s) Time step
%   Q          - 过程噪声协方差矩阵 (6x6) Process noise covariance matrix
%
% Outputs:
%   x_pred     - 预测状态 [q; b] (7x1) Predicted state
%   P_pred     - 预测协方差 (6x6) Predicted covariance

    % 解包状态 Unpack state
    q_hat = x_est(1:4);
    b_hat = x_est(5:7);

    % 1. 估计真实角速度 Estimate True Angular Velocity
    % omega_hat = omega_meas - b_hat
    w_hat = omega_meas - b_hat;
    w_norm = norm(w_hat);

    % 2. 四元数传播 (闭式解) Propagate Quaternion (Closed Form Solution)
    % dq/dt = 0.5 * Omega(w) * q
    % For constant w, q(t+dt) = exp(0.5*Omega*dt) * q(t)
    
    if w_norm > 1e-8
        c = cos(0.5 * w_norm * dt);
        s = sin(0.5 * w_norm * dt);
        Om = utils_Omega(w_hat);
        
        % 状态转移矩阵 (四元数部分) Transition matrix for quaternion (4x4)
        Phi_q = c * eye(4) + (s / w_norm) * Om;
    else
        % 小角度近似 Small angle approximation
        Om = utils_Omega(w_hat);
        Phi_q = eye(4) + 0.5 * dt * Om;
    end
    
    q_pred = Phi_q * q_hat;
    q_pred = q_pred / norm(q_pred); % 强制归一化 Enforce normalization

    % 3. 漂移传播 Propagate Bias
    % b_dot = 0 (random walk mean)
    b_pred = b_hat;

    % 重组状态 Reassemble state
    x_pred = [q_pred; b_pred];

    % 4. 误差状态转移矩阵 (6x6) State Transition Matrix for Error State
    % Delta x = [Delta theta; Delta b]
    % Phi = [ Theta   Psi ]
    %       [ 0       I   ]
    
    % Theta: 姿态误差转移矩阵
    % Theta approx = I - [w x]*dt
    
    wx = [ 0, -w_hat(3), w_hat(2);
           w_hat(3), 0, -w_hat(1);
          -w_hat(2), w_hat(1), 0 ];
          
    % 一阶近似 First order approximation
    Theta = eye(3) - wx * dt;
    
    % Psi: 偏差误差对姿态误差的耦合
    % Psi = -I * dt (approx)
    Psi = -eye(3) * dt;
    
    Phi = [ Theta, Psi;
            zeros(3), eye(3) ];

    % 5. 协方差传播 Propagate Covariance
    P_pred = Phi * P_est * Phi' + Q;

end
