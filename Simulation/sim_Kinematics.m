function [time, q_true, w_true] = sim_Kinematics(T_total, dt)
% SIM_KINEMATICS 模拟真实姿态动力学 Simulates the ground truth attitude dynamics
%
% 生成航天器姿态轨迹 Generates a trajectory for the spacecraft attitude.
%
% Inputs:
%   T_total - 总仿真时间 (秒) Total simulation time (seconds)
%   dt      - 采样间隔 (秒) Samping internal (seconds)
%
% Outputs:
%   time   - 时间向量 Time vector
%   q_true - 真实四元数 True quaternions (4 x N)
%   w_true - 真实角速度 True angular velocity (3 x N)

    % 时间设置 Time setup
    time = 0:dt:T_total;
    N = length(time);
    
    % 初始化数组 Initialize arrays
    q_true = zeros(4, N);
    w_true = zeros(3, N);
    
    % 初始条件 Initial conditions
    % 从单位四元数开始 (与惯性系对齐) Start with identity quaternion (Aligned with inertial frame)
    q_curr = [0; 0; 0; 1];
    q_true(:, 1) = q_curr;
    
    % 定义角速度曲线 (真值) Define Angular Velocity Profile (Ground Truth)
    % 示例：各轴的正弦运动 Example: Sinusoidal motion on all axes
    for k = 1:N
        t = time(k);
        % 动态运动： A dynamic motion:
        w_curr = [ 0.1 * sin(0.1*t);
                   0.05 * cos(0.2*t);
                   0.02 * sin(0.05*t) + 0.02 ];
        w_true(:, k) = w_curr;
    end
    
    % 积分循环 (RK4) Integration Loop (RK4)
    for k = 1:N-1
        w_now = w_true(:, k);
        w_next = w_true(:, k+1); % 简单近似，或全步使用 w_now (Simple approximation, or use w_now for whole step)
        
        q_now = q_true(:, k);
        
        % 由于 w 变化缓慢，使用简单积分 Use simple integration since w is changing slowly
        % 如果需要，可以使用带插值 w 的 RK4 (Or use RK4 with interpolated w if needed.)
        % 这里假设 w 在 dt 步长内恒定 (零阶保持) Here we assume w is constant over dt for the step (Zero Order Hold)
        
        k1 = 0.5 * utils_Omega(w_now) * q_now;
        
        q_temp2 = q_now + 0.5 * dt * k1;
        k2 = 0.5 * utils_Omega(w_now) * q_temp2; % 使用 w_now 近似 Using w_now approximation
        
        q_temp3 = q_now + 0.5 * dt * k2;
        k3 = 0.5 * utils_Omega(w_now) * q_temp3;
        
        q_temp4 = q_now + dt * k3;
        k4 = 0.5 * utils_Omega(w_now) * q_temp4;
        
        q_next = q_now + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
        
        % 归一化四元数 Normalize Quaternion
        q_next = q_next / norm(q_next);
        
        q_true(:, k+1) = q_next;
    end

end
