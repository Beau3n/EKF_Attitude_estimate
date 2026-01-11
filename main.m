% main.m
% Attitude Determination System EKF Simulation
% 姿态确定系统 EKF 仿真主程序
%
% This script sets up the simulation environment, generates ground truth
% trajectories, simulates sensor measurements (Gyro & Star Tracker),
% runs the Extended Kalman Filter (EKF), and analyzes the performance.

clc; clear; close all;

% Add paths
addpath(genpath('.'));

%% 1. Simulation Configuration (仿真配置)

% Time settings
dt = 0.1;               % Sampling interval (s) (10 Hz)
T_total = 600;          % Total simulation time (s)
time = 0:dt:T_total;
N = length(time);

% Sensor Parameters
% 1. Gyroscope (陀螺仪)
params_gyro.gyro_noise_std = deg2rad(0.01);    % Noise (ARW) [rad/s]
params_gyro.gyro_bias_std  = deg2rad(1e-5);    % Bias Instability (RRW) [rad/s/sqrt(s)] (Simplified usage)
params_gyro.b_init         = deg2rad([0.1; -0.2; 0.1]); % True Initial Bias [rad/s]

% 2. Star Tracker (星敏感器)
params_st.st_noise_std     = deg2rad(0.005);   % Measurement Noise (18 arcsec) [rad]

%% 2. Ground Truth Generation (真值生成)

fprintf('Generating True Trajectory...\n');
[time, q_true, w_true] = sim_Kinematics(T_total, dt);

% Calculate True Euler Angles (ZYX) for reference
euler_true = zeros(3, N);
for k = 1:N
    % Need a quat2eul function or implement here.
    % Implementing simple ZYX conversion
    q = q_true(:, k);
    % q = [q1, q2, q3, q4] (scalar last format used in this project?)
    % Let's verify utils... Predict.m uses q(1:4).
    % Assuming q = [v; s] format based on Predict.m: q_pred = Phi_q * q_hat.
    % And utils_Omega: q4 is scalar part. YES.
    
    % Standard ZYX Euler (Roll, Pitch, Yaw)
    % R = utils_AttitudeMatrix(q);
    % [yaw, pitch, roll] can be extracted from R
    % For simplicity in error analysis, we usually look at error angles, not absolute euler differences.
end

%% 3. Sensor Simulation (传感器仿真)

fprintf('Simulating Sensor Measurements...\n');

% 3.1 Gyroscope
[gyro_meas, b_true] = Model_Gyro_Sensor(w_true, dt, params_gyro);

% 3.2 Star Tracker
% Define Reference Stars in Inertial Frame (Unit vectors)
% Using 2 stars for full observability
v_ref = [1, 0, 0; 0, 1, 0; 0, 0, 1]'; % 3 orthogonal stars to ensure observability
% Or just 2 stars separated by some angle
v_ref = [1, 0, 0; 
         0, 0, 1]'; 

[star_meas] = Model_StarTracker(q_true, v_ref, params_st);

%% 4. EKF Initialization (滤波器初始化)

fprintf('Initializing EKF...\n');

% Initial State Guess
% Add large initial error to test convergence
delta_angle_init = deg2rad(10); % 10 degrees initial error
axis_init = [1; 1; 1] / sqrt(3);
dq_init = [axis_init * sin(delta_angle_init/2); cos(delta_angle_init/2)];

% q_init = q_true(:,1) * dq_init (quaternion multiplication)
% Note: q_init = q_true x dq_init
q_init_true = q_true(:, 1);
% q_init = q_multiply(q_init_true, dq_init); % We don't have a lib function yet
% Let's just assign identity or a fixed offset for simplicity if multiply isn't handy
% Using small angle approximation: q_new approx q_true + 0.5 * Xi * dtheta
q_est_init = [0; 0; 0; 1]; % Initialize with Identity (Cold start)

b_est_init = [0; 0; 0];    % Initialize bias as 0

x_init = [q_est_init; b_est_init];

% Initial Covariance P (6x6)
% High uncertainty for attitude, low for bias
P_init = zeros(6);
P_init(1:3, 1:3) = (deg2rad(20))^2 * eye(3); % Attitude uncertainty
P_init(4:6, 4:6) = (deg2rad(0.1))^2 * eye(3); % Bias uncertainty

% Process Noise Covariance Q (6x6)
% Tuned based on sensor specs
sig_v = params_gyro.gyro_noise_std;
sig_u = params_gyro.gyro_bias_std;

% Continuous Q spectral density
% Q_cont = [sig_v^2 * I,  0;
%           0,            sig_u^2 * I]
% Discrete approximation: Q_k = Q_cont * dt
Q = zeros(6);
scale_factor = 1.0; % Tuning factor
Q(1:3, 1:3) = scale_factor * (sig_v^2 * dt + 1/3 * sig_u^2 * dt^3) * eye(3);
Q(4:6, 4:6) = scale_factor * (sig_u^2 * dt) * eye(3);

% Measurement Noise Covariance R (3x3 per vector)
R = (params_st.st_noise_std^2) * eye(3);

%% 5. Run EKF (运行滤波)

fprintf('Running EKF...\n');
[x_hist, P_hist] = EKF_update(time, gyro_meas, star_meas, v_ref, x_init, P_init, Q, R);

%% 6. Analysis & Plotting (结果分析与绘图)

fprintf('Analyzing Results...\n');

% 6.1 Compute Attitude Error (delta theta)
% Error quaternion: q_err = q_true^(-1) * q_est
% q_inv = [-v; s]
err_euler_deg = zeros(3, N);
err_sigma_deg = zeros(3, N);

for k = 1:N
    q_t = q_true(:, k);
    q_e = x_hist(1:4, k);
    
    % Inverse of true quaternion
    q_t_inv = [-q_t(1:3); q_t(4)];
    
    % Quaternion multiplication: q_err = q_t_inv * q_e
    % Formula:
    % [s1*v2 + s2*v1 + v1 x v2]
    % [s1*s2 - v1.v2]
    
    xyz_t = q_t_inv(1:3); s_t = q_t_inv(4);
    xyz_e = q_e(1:3);     s_e = q_e(4);
    
    xyz_err = s_t * xyz_e + s_e * xyz_t + cross(xyz_t, xyz_e);
    s_err   = s_t * s_e - dot(xyz_t, xyz_e);
    
    % Convert to small angle errors (in body frame)
    % delta_theta = 2 * vector_part
    % If scalar part is negative, invert to keep closest path
    if s_err < 0
        xyz_err = -xyz_err;
        s_err = -s_err;
    end
    
    d_theta = 2 * xyz_err;
    err_euler_deg(:, k) = rad2deg(d_theta);
    
    % Extract 3-sigma bounds
    P_k = P_hist(:, :, k);
    diag_P = diag(P_k);
    std_att = sqrt(diag_P(1:3));
    err_sigma_deg(:, k) = 3 * rad2deg(std_att);
end

% 6.2 Compute Bias Error
b_est_hist = x_hist(5:7, :);
b_err = b_est_hist - b_true;
b_err_deg = rad2deg(b_err);
b_sigma_deg = zeros(3, N);

for k = 1:N
    P_k = P_hist(:, :, k);
    diag_P = diag(P_k);
    std_b = sqrt(diag_P(4:6));
    b_sigma_deg(:, k) = 3 * rad2deg(std_b);
end

% 6.3 Performance Metrics
% Convergence time: time to stay within bounds (e.g., 0.1 deg)
% Steady state RMSE (Last 30% of data)
idx_ss = floor(0.7 * N):N;
rmse_att = sqrt(mean(err_euler_deg(:, idx_ss).^2, 2));
rmse_bias = sqrt(mean(b_err_deg(:, idx_ss).^2, 2));

fprintf('\n=== Performance Metrics ===\n');
fprintf('Simulation Time: %.1f s\n', T_total);
fprintf('Steady State RMSE (Attitude) [deg]:\n  X: %.4f, Y: %.4f, Z: %.4f\n', rmse_att);
fprintf('Steady State RMSE (Bias) [deg/s]:\n  X: %.4e, Y: %.4e, Z: %.4e\n', rmse_bias);

%% 7. Plotting

% Plot Attitude Errors
figure('Name', 'Attitude Estimation Error', 'Color', 'w');
titles = {'Roll Error (X)', 'Pitch Error (Y)', 'Yaw Error (Z)'};
for i = 1:3
    subplot(3, 1, i); hold on; grid on;
    % Plot error
    plot(time, err_euler_deg(i, :), 'b', 'LineWidth', 1.0);
    % Plot 3-sigma bounds
    plot(time, err_sigma_deg(i, :), 'r--', 'LineWidth', 1.0);
    plot(time, -err_sigma_deg(i, :), 'r--', 'LineWidth', 1.0);
    
    ylabel('Error [deg]');
    title(titles{i});
    if i == 1
        legend('Error', '3\sigma Bound');
    end
end
xlabel('Time [s]');

% Plot Bias Estimation
figure('Name', 'Gyro Bias Estimation', 'Color', 'w');
titles_b = {'Bias X', 'Bias Y', 'Bias Z'};
for i = 1:3
    subplot(3, 1, i); hold on; grid on;
    plot(time, rad2deg(b_est_hist(i, :)), 'b', 'LineWidth', 1.2);
    plot(time, rad2deg(b_true(i, :)), 'g--', 'LineWidth', 1.2);
    ylabel('Bias [deg/s]');
    legend('Estimated', 'True');
    title(titles_b{i});
end
xlabel('Time [s]');

% Plot Bias Error
figure('Name', 'Bias Estimation Error', 'Color', 'w');
for i = 1:3
    subplot(3, 1, i); hold on; grid on;
    plot(time, b_err_deg(i, :), 'b', 'LineWidth', 1.0);
    plot(time, b_sigma_deg(i, :), 'r--', 'LineWidth', 1.0);
    plot(time, -b_sigma_deg(i, :), 'r--', 'LineWidth', 1.0);
    ylabel('Error [deg/s]');
    title(['Bias Error ' titles_b{i}]);
end
xlabel('Time [s]');
