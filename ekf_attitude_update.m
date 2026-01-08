function x_new = ekf_attitude_update(x_est, omega_meas, q_star, dt, Q, R)
% ekf_attitude_update 基于星敏感器/陀螺组合的EKF姿态确定算法
%
% 采用乘性扩展卡尔曼滤波 (MEKF) 方法
% 状态向量: 误差状态 δx = [δθ(3x1); δβ(3x1)]
%   δθ: 姿态误差角 (机体系)
%   δβ: 陀螺仪零偏误差
%
% 输入:
%   x_est : 当前状态估计结构体
%       .q    : 姿态四元数 (4x1) [q0; q1; q2; q3]
%       .beta : 陀螺仪零偏估计 (3x1) [rad/s]
%       .P    : 误差协方差矩阵 (6x6)
%   omega_meas : 陀螺仪测量角速度 (3x1) [rad/s]
%   q_star     : 星敏感器测量四元数 (4x1)
%   dt         : 采样周期 [s]
%   Q          : 过程噪声协方差矩阵 (6x6)
%   R          : 测量噪声协方差矩阵 (3x3)
%
% 输出:
%   x_new : 更新后的状态估计结构体

    % 提取当前状态
    q_k = x_est.q;
    beta_k = x_est.beta;
    P_k = x_est.P;
    
    %% ========== 1. 时间更新 (预测) ==========
    
    % 1.1 修正角速度 (去除零偏估计)
    omega_hat = omega_meas - beta_k;
    
    % 1.2 四元数传播 (姿态外推)
    % 使用精确的四元数运动学方程离散化
    omega_norm = norm(omega_hat);
    
    if omega_norm > 1e-12
        % 旋转增量四元数
        theta = omega_norm * dt;
        axis_n = omega_hat / omega_norm;
        dq = [cos(theta/2); axis_n * sin(theta/2)];
    else
        dq = [1; 0; 0; 0];
    end
    
    % 四元数更新: q(k+1) = q(k) ⊗ dq
    q_pred = quat_multiply(q_k, dq);
    q_pred = q_pred / norm(q_pred);  % 归一化
    
    % 1.3 零偏预测 (随机游走模型，预测值不变)
    beta_pred = beta_k;
    
    % 1.4 状态转移矩阵 Φ
    % 连续时间状态方程:
    %   δθ_dot = -[ω×]δθ - δβ + η_v
    %   δβ_dot = η_u
    % 离散化: Φ ≈ I + F*dt
    
    Omega_cross = skew_symmetric(omega_hat);
    
    F = [-Omega_cross, -eye(3);
         zeros(3),     zeros(3)];
    
    Phi = eye(6) + F * dt;
    
    % 1.5 协方差预测
    P_pred = Phi * P_k * Phi' + Q;
    
    %% ========== 2. 测量更新 (校正) ==========
    
    % 2.1 计算测量残差 (误差四元数)
    % δq = q_star ⊗ q_pred^(-1)
    q_pred_inv = [q_pred(1); -q_pred(2:4)];  % 四元数共轭即为逆
    dq_meas = quat_multiply(q_star, q_pred_inv);
    
    % 确保标量部分为正 (避免符号问题)
    if dq_meas(1) < 0
        dq_meas = -dq_meas;
    end
    
    % 2.2 提取小角度误差作为观测残差
    % δθ ≈ 2 * [δq1; δq2; δq3]
    z = 2 * dq_meas(2:4);
    
    % 2.3 观测矩阵 H
    % z = H * δx + v
    % 对于四元数直接观测: H = [I_3x3, 0_3x3]
    H = [eye(3), zeros(3)];
    
    % 2.4 卡尔曼增益
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    
    % 2.5 状态修正
    dx = K * z;
    d_theta = dx(1:3);
    d_beta = dx(4:6);
    
    % 2.6 四元数更新 (乘性修正)
    % q_new = q_pred ⊗ δq(d_theta)
    dq_corr = [1; 0.5 * d_theta];
    dq_corr = dq_corr / norm(dq_corr);
    q_new = quat_multiply(q_pred, dq_corr);
    q_new = q_new / norm(q_new);
    
    % 2.7 零偏更新
    beta_new = beta_pred + d_beta;
    
    % 2.8 协方差更新 (Joseph形式，保证正定性)
    I_KH = eye(6) - K * H;
    P_new = I_KH * P_pred * I_KH' + K * R * K';
    
    %% ========== 3. 输出 ==========
    x_new.q = q_new;
    x_new.beta = beta_new;
    x_new.P = P_new;

end

%% ==================== 辅助函数 ====================

function q = quat_multiply(q1, q2)
% 四元数乘法 q = q1 ⊗ q2
% q = [q0; q1; q2; q3], q0为标量
    s1 = q1(1); v1 = q1(2:4);
    s2 = q2(1); v2 = q2(2:4);
    
    s = s1*s2 - dot(v1, v2);
    v = s1*v2 + s2*v1 + cross(v1, v2);
    q = [s; v(:)];
end

function S = skew_symmetric(v)
% 反对称矩阵 [v×]
    S = [  0,   -v(3),  v(2);
         v(3),   0,   -v(1);
        -v(2),  v(1),   0  ];
end
