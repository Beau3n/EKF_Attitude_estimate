function [omega_meas, bias_true] = get_gyro_mode(omega_true, bias_true, dt, params)
% get_gyro_mode 陀螺仪测量模型
% 
% 输入:
%   omega_true : 真实角速度 (3x1) [rad/s]
%   bias_true  : 当前真实零偏 (3x1) [rad/s]
%   dt         : 采样间隔 [s]
%   params     : 陀螺仪参数结构体
%       .sigma_v : 角度随机游走系数 (ARW) [rad/s/sqrt(Hz)]
%       .sigma_u : 零偏不稳定性系数 (RRW) [rad/s^2/sqrt(Hz)]
%
% 输出:
%   omega_meas : 陀螺仪测量输出 (3x1) [rad/s]
%   bias_true  : 更新后的真实零偏 (3x1) [rad/s]

    % 提取参数
    sigma_v = params.sigma_v;  % 角度随机游走
    sigma_u = params.sigma_u;  % 零偏随机游走
    
    % 零偏演化 (随机游走模型)
    % bias(k+1) = bias(k) + sigma_u * sqrt(dt) * w_u
    bias_true = bias_true + sigma_u * sqrt(dt) * randn(3, 1);
    
    % 测量噪声 (白噪声)
    % noise = sigma_v / sqrt(dt) * w_v
    noise = sigma_v / sqrt(dt) * randn(3, 1);
    
    % 陀螺仪测量模型: omega_meas = omega_true + bias + noise
    omega_meas = omega_true + bias_true + noise;

end
