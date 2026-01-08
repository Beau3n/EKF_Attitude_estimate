function q_meas = star_tracker(q_true, params)
% star_tracker 星敏感器测量模型
% 
% 星敏感器直接输出姿态四元数，带有测量噪声
%
% 输入:
%   q_true : 真实姿态四元数 (4x1) [q0; q1; q2; q3], q0为标量
%   params : 星敏感器参数结构体
%       .sigma_st : 星敏感器测量噪声标准差 [rad] (各轴)
%
% 输出:
%   q_meas : 星敏感器测量四元数 (4x1)

    sigma_st = params.sigma_st;
    
    % 生成小角度误差 (三轴)
    delta_theta = sigma_st .* randn(3, 1);
    
    % 将小角度误差转换为误差四元数
    % delta_q ≈ [1; 0.5*delta_theta] (小角度近似)
    delta_q = [1; 0.5 * delta_theta];
    delta_q = delta_q / norm(delta_q);
    
    % 测量四元数 = 真实四元数 ⊗ 误差四元数
    q_meas = quat_multiply(q_true, delta_q);
    q_meas = q_meas / norm(q_meas);

end

%% 辅助函数: 四元数乘法
function q = quat_multiply(q1, q2)
    % q = q1 ⊗ q2
    % 四元数定义: q = [q0; q1; q2; q3], q0为标量部分
    s1 = q1(1); v1 = q1(2:4);
    s2 = q2(1); v2 = q2(2:4);
    
    s = s1*s2 - dot(v1, v2);
    v = s1*v2 + s2*v1 + cross(v1, v2);
    q = [s; v(:)];
end
