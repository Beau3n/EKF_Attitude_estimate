% Test_Prediction.m
% 测试 EKF 预测步骤 (纯陀螺积分) vs 真值
% Test EKF Prediction step (Pure Gyro Integration) vs Ground Truth

clc; clear; close all;
addpath('Utils');
addpath('Simulation');
addpath('EKF');

%% 1. 配置 Configuration
dt = 0.01;       % 高频采样 High sampling rate
T_total = 600;    % 60秒仿真 60s simulation

% 传感器参数 (包含偏差)
params.gyro_noise_std = 1e-5; % 噪声 Noise
params.gyro_bias_std  = 1e-6; 
params.b_init = [0.005; -0.005; 0.002]; % 显著的初始偏差 Significant initial bias

%% 2. 生成数据 Generate Data
fprintf('生成仿真数据 Generating Data...\n');
[time, q_true, w_true] = sim_Kinematics(T_total, dt);
[gyro_meas, b_true] = Model_Gyro_Sensor(w_true, dt, params);

%% 3. 运行预测 Testing Prediction
fprintf('运行预测测试 Running Prediction Test...\n');
N = length(time);

% 初始化估计 Initialization
% 假设初始姿态已知 Perfect initial attitude
x_est = [q_true(:, 1); 0; 0; 0]; % 初始 Bias 设为 0 (未校准) Initial Bias set to 0 (Uncalibrated)
P_est = eye(6) * 1e-6;
Q = eye(6) * 1e-6; % 随意设置，此处主要看状态传播 Set arbitrarily, mainly looking at state propagation

q_pred_hist = zeros(4, N);
angle_error = zeros(1, N);

for k = 1:N
    % 运行预测 Run Predict
    if k > 1
        x_est = EKF_prediction(x_est, P_est, gyro_meas(:, k-1), dt, Q);
    end
    
    % 记录 Record
    q_pred = x_est(1:4);
    q_pred_hist(:, k) = q_pred;
    
    % 计算与真值的误差角度 Calculate error angle vs Truth
    % error = 2 * acos( |q_pred . q_true| )
    dot_prod = abs(q_pred' * q_true(:, k));
    if dot_prod > 1, dot_prod = 1; end
    angle_error(k) = 2 * acos(dot_prod);
end

%% 4. 绘图结果 Visualization
figure('Name', 'Prediction Drift Test');

% 绘制欧拉角或者四元数分量 Comparison of Quaternion
subplot(2,1,1);
plot(time, q_true(1,:), 'k-', 'LineWidth', 1.5); hold on;
plot(time, q_pred_hist(1,:), 'r--');
title('Quaternion q1: True vs Predicted (Gyro Integration)');
legend('True', 'Predicted');
grid on;

% 绘制角度误差 Plot Angle Error
subplot(2,1,2);
plot(time, rad2deg(angle_error), 'r', 'LineWidth', 1.5);
title(['Attitude Error Growth (Drift due to Uncompensated Bias: ' mat2str(params.b_init') ')']);
xlabel('Time (s)');
ylabel('Error (deg)');
grid on;

fprintf('测试完成。由于未校正陀螺零偏，预期误差会随时间发散。\nTest Complete. Error is expected to diverge due to uncompensated gyro bias.\n');
