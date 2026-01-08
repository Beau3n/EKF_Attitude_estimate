%% 基于EKF的星敏感器/陀螺组合姿态确定仿真
% 参考论文: 基于EKF的航天器姿态确定算法及精度分析
% 
% 传感器配置:
%   - 陀螺仪: 提供角速度测量 (高采样率)
%   - 星敏感器: 提供姿态四元数测量 (低采样率)
%
% EKF状态向量:
%   - 姿态四元数 q (4x1)
%   - 陀螺仪零偏 β (3x1)

clear; clc; close all;

%% ==================== 1. 仿真参数设置 ====================

% 时间参数
T_total = 600;          % 总仿真时间 [s]
dt_gyro = 0.01;         % 陀螺仪采样周期 [s] (100 Hz)
dt_star = 0.5;          % 星敏感器采样周期 [s] (2 Hz)
time = 0:dt_gyro:T_total;
N = length(time);

% 陀螺仪参数 (典型光纤陀螺)
gyro_params.sigma_v = 1e-4;    % 角度随机游走 [rad/s/sqrt(Hz)]
gyro_params.sigma_u = 1e-6;    % 零偏不稳定性 [rad/s^2/sqrt(Hz)]

% 星敏感器参数
star_params.sigma_st = [5; 5; 30] * (1/3600) * (pi/180);  % 测量噪声 [rad]
% 横滚/俯仰: 5角秒, 偏航: 30角秒 (典型星敏精度)

% 真实初始零偏
bias_true_init = [0.5; -0.3; 0.2] * (pi/180) / 3600;  % [rad/s] (约0.5°/h)

%% ==================== 2. 真值轨迹生成 ====================

% 真实初始姿态 (单位四元数)
q_true = [1; 0; 0; 0];

% 航天器角速度 (模拟三轴稳定卫星，含周期性扰动)
omega_true_func = @(t) [0.001*sin(0.01*t); 
                        0.002*cos(0.01*t); 
                        0.0005];  % [rad/s]

% 预分配存储
q_true_history = zeros(4, N);
omega_true_history = zeros(3, N);
bias_true_history = zeros(3, N);

% 真实零偏 (初始化)
bias_true = bias_true_init;

%% ==================== 3. EKF 初始化 ====================

% 初始姿态估计 (引入初始误差)
init_err_angle = 5 * (pi/180);  % 初始姿态误差 5°
init_err_axis = [1; 1; 1] / sqrt(3);
dq_init_err = [cos(init_err_angle/2); init_err_axis * sin(init_err_angle/2)];
x_est.q = quat_multiply(q_true, dq_init_err);
x_est.q = x_est.q / norm(x_est.q);

% 初始零偏估计 (假设未知，设为0)
x_est.beta = [0; 0; 0];

% 初始协方差矩阵
P0_att = (10 * pi/180)^2 * eye(3);     % 姿态初始不确定度 10°
P0_bias = (1 * pi/180/3600)^2 * eye(3); % 零偏初始不确定度 1°/h
x_est.P = blkdiag(P0_att, P0_bias);

% 过程噪声协方差矩阵 Q
Q_att = (gyro_params.sigma_v)^2 * dt_gyro * eye(3);
Q_bias = (gyro_params.sigma_u)^2 * dt_gyro * eye(3);
Q = blkdiag(Q_att, Q_bias);

% 测量噪声协方差矩阵 R (星敏感器)
R = diag(star_params.sigma_st.^2);

%% ==================== 4. 数据存储 ====================

q_est_history = zeros(4, N);
bias_est_history = zeros(3, N);
att_err_history = zeros(3, N);  % 姿态误差 [rad]
bias_err_history = zeros(3, N); % 零偏误差 [rad/s]
P_diag_history = zeros(6, N);   % 协方差对角线

%% ==================== 5. 主仿真循环 ====================

fprintf('开始仿真...\n');
star_update_counter = 0;

for k = 1:N
    t = time(k);
    
    %% 5.1 真值生成
    omega_true = omega_true_func(t);
    
    % 真实姿态传播
    if k > 1
        [q_true, bias_true] = propagate_true_state(q_true, omega_true, ...
                                                    bias_true, dt_gyro, gyro_params);
    end
    
    % 存储真值
    q_true_history(:, k) = q_true;
    omega_true_history(:, k) = omega_true;
    bias_true_history(:, k) = bias_true;
    
    %% 5.2 传感器测量生成
    % 陀螺仪测量 (每步都有)
    [omega_meas, bias_true] = get_gyro_mode(omega_true, bias_true, dt_gyro, gyro_params);
    
    % 星敏感器测量 (周期性)
    star_update_counter = star_update_counter + dt_gyro;
    has_star_update = (star_update_counter >= dt_star);
    
    if has_star_update
        q_star = star_tracker(q_true, star_params);
        star_update_counter = 0;
    else
        q_star = [];
    end
    
    %% 5.3 EKF 时间更新 (每个陀螺采样周期)
    x_est = ekf_time_update(x_est, omega_meas, dt_gyro, Q);
    
    %% 5.4 EKF 测量更新 (有星敏数据时)
    if has_star_update && ~isempty(q_star)
        x_est = ekf_measurement_update(x_est, q_star, R);
    end
    
    %% 5.5 记录估计结果
    q_est_history(:, k) = x_est.q;
    bias_est_history(:, k) = x_est.beta;
    P_diag_history(:, k) = diag(x_est.P);
    
    % 计算误差
    [att_err, ~] = compute_attitude_error(x_est.q, q_true);
    att_err_history(:, k) = att_err;
    bias_err_history(:, k) = x_est.beta - bias_true;
    
end

fprintf('仿真完成!\n');

%% ==================== 6. 结果可视化 ====================

% 转换单位
att_err_deg = att_err_history * (180/pi);           % 转为度
att_err_arcsec = att_err_history * (180/pi) * 3600; % 转为角秒
bias_err_deg_h = bias_err_history * (180/pi) * 3600; % 转为°/h
P_att_3sigma = 3 * sqrt(P_diag_history(1:3, :)) * (180/pi) * 3600; % 3σ边界

%% 图1: 姿态估计误差
figure('Name', '姿态估计误差', 'Position', [100, 100, 900, 600]);

subplot(3,1,1);
plot(time, att_err_arcsec(1,:), 'b', 'LineWidth', 0.5);
hold on;
plot(time, P_att_3sigma(1,:), 'r--', time, -P_att_3sigma(1,:), 'r--', 'LineWidth', 1);
ylabel('Roll (角秒)');
title('姿态估计误差 (3\sigma 边界)');
legend('估计误差', '3\sigma边界');
grid on;

subplot(3,1,2);
plot(time, att_err_arcsec(2,:), 'b', 'LineWidth', 0.5);
hold on;
plot(time, P_att_3sigma(2,:), 'r--', time, -P_att_3sigma(2,:), 'r--', 'LineWidth', 1);
ylabel('Pitch (角秒)');
grid on;

subplot(3,1,3);
plot(time, att_err_arcsec(3,:), 'b', 'LineWidth', 0.5);
hold on;
plot(time, P_att_3sigma(3,:), 'r--', time, -P_att_3sigma(3,:), 'r--', 'LineWidth', 1);
ylabel('Yaw (角秒)');
xlabel('时间 (s)');
grid on;

%% 图2: 陀螺零偏估计
figure('Name', '陀螺零偏估计', 'Position', [150, 150, 900, 600]);

subplot(3,1,1);
plot(time, bias_true_history(1,:)*(180/pi)*3600, 'k--', 'LineWidth', 1.5);
hold on;
plot(time, bias_est_history(1,:)*(180/pi)*3600, 'b', 'LineWidth', 1);
ylabel('X轴 (°/h)');
title('陀螺零偏估计');
legend('真值', '估计值');
grid on;

subplot(3,1,2);
plot(time, bias_true_history(2,:)*(180/pi)*3600, 'k--', 'LineWidth', 1.5);
hold on;
plot(time, bias_est_history(2,:)*(180/pi)*3600, 'b', 'LineWidth', 1);
ylabel('Y轴 (°/h)');
grid on;

subplot(3,1,3);
plot(time, bias_true_history(3,:)*(180/pi)*3600, 'k--', 'LineWidth', 1.5);
hold on;
plot(time, bias_est_history(3,:)*(180/pi)*3600, 'b', 'LineWidth', 1);
ylabel('Z轴 (°/h)');
xlabel('时间 (s)');
grid on;

%% 图3: 零偏估计误差
figure('Name', '零偏估计误差', 'Position', [200, 200, 900, 400]);
plot(time, bias_err_deg_h', 'LineWidth', 1);
ylabel('误差 (°/h)');
xlabel('时间 (s)');
title('陀螺零偏估计误差');
legend('X轴', 'Y轴', 'Z轴');
grid on;

%% 统计结果
fprintf('\n========== 仿真统计结果 ==========\n');
fprintf('姿态估计精度 (稳态RMS, 后50%%数据):\n');
idx_steady = round(N/2):N;
rms_att = rms(att_err_arcsec(:, idx_steady), 2);
fprintf('  Roll:  %.2f 角秒\n', rms_att(1));
fprintf('  Pitch: %.2f 角秒\n', rms_att(2));
fprintf('  Yaw:   %.2f 角秒\n', rms_att(3));

rms_bias = rms(bias_err_deg_h(:, idx_steady), 2);
fprintf('\n零偏估计精度 (稳态RMS):\n');
fprintf('  X: %.4f °/h\n', rms_bias(1));
fprintf('  Y: %.4f °/h\n', rms_bias(2));
fprintf('  Z: %.4f °/h\n', rms_bias(3));

%% ==================== 辅助函数 ====================

function [q_new, bias_new] = propagate_true_state(q, omega, bias, dt, params)
% 真值状态传播
    omega_norm = norm(omega);
    if omega_norm > 1e-12
        theta = omega_norm * dt;
        n = omega / omega_norm;
        dq = [cos(theta/2); n*sin(theta/2)];
    else
        dq = [1; 0; 0; 0];
    end
    q_new = quat_multiply(q, dq);
    q_new = q_new / norm(q_new);
    
    % 零偏随机游走
    bias_new = bias + params.sigma_u * sqrt(dt) * randn(3,1);
end

function x_new = ekf_time_update(x_est, omega_meas, dt, Q)
% EKF 时间更新
    q_k = x_est.q;
    beta_k = x_est.beta;
    P_k = x_est.P;
    
    % 角速度修正
    omega_hat = omega_meas - beta_k;
    
    % 四元数传播
    omega_norm = norm(omega_hat);
    if omega_norm > 1e-12
        theta = omega_norm * dt;
        n = omega_hat / omega_norm;
        dq = [cos(theta/2); n*sin(theta/2)];
    else
        dq = [1; 0; 0; 0];
    end
    q_pred = quat_multiply(q_k, dq);
    q_pred = q_pred / norm(q_pred);
    
    % 状态转移矩阵
    Omega_cross = skew_symmetric(omega_hat);
    F = [-Omega_cross, -eye(3); zeros(3), zeros(3)];
    Phi = eye(6) + F * dt;
    
    % 协方差预测
    P_pred = Phi * P_k * Phi' + Q;
    
    x_new.q = q_pred;
    x_new.beta = beta_k;
    x_new.P = P_pred;
end

function x_new = ekf_measurement_update(x_est, q_star, R)
% EKF 测量更新
    q_pred = x_est.q;
    beta_pred = x_est.beta;
    P_pred = x_est.P;
    
    % 误差四元数
    q_pred_inv = [q_pred(1); -q_pred(2:4)];
    dq_meas = quat_multiply(q_star, q_pred_inv);
    if dq_meas(1) < 0
        dq_meas = -dq_meas;
    end
    
    % 观测残差
    z = 2 * dq_meas(2:4);
    
    % 观测矩阵
    H = [eye(3), zeros(3)];
    
    % 卡尔曼增益
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    
    % 状态修正
    dx = K * z;
    d_theta = dx(1:3);
    d_beta = dx(4:6);
    
    % 四元数更新
    dq_corr = [1; 0.5*d_theta];
    dq_corr = dq_corr / norm(dq_corr);
    q_new = quat_multiply(q_pred, dq_corr);
    q_new = q_new / norm(q_new);
    
    % 零偏更新
    beta_new = beta_pred + d_beta;
    
    % 协方差更新 (Joseph形式)
    I_KH = eye(6) - K*H;
    P_new = I_KH * P_pred * I_KH' + K * R * K';
    
    x_new.q = q_new;
    x_new.beta = beta_new;
    x_new.P = P_new;
end

function [att_err, dq] = compute_attitude_error(q_est, q_true)
% 计算姿态误差角
    q_est_inv = [q_est(1); -q_est(2:4)];
    dq = quat_multiply(q_true, q_est_inv);
    if dq(1) < 0
        dq = -dq;
    end
    att_err = 2 * dq(2:4);
end

function q = quat_multiply(q1, q2)
% 四元数乘法
    s1 = q1(1); v1 = q1(2:4);
    s2 = q2(1); v2 = q2(2:4);
    s = s1*s2 - dot(v1, v2);
    v = s1*v2 + s2*v1 + cross(v1, v2);
    q = [s; v(:)];
end

function S = skew_symmetric(v)
% 反对称矩阵
    S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end
