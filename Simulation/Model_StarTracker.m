function star_meas = Model_StarTracker(q_true, v_ref, params)
% MODEL_STARTRACKER Simulates Star Tracker Vector Measurements
%
% Model: z = A(q)*r + noise
%
% Inputs:
%   q_true - True quaternion history (4 x N)
%   v_ref  - Inertial reference vectors (3 x M), M is number of stars
%   params - Struct with fields:
%            params.st_noise_std (Measurement noise std dev in rad or unitless)
%
% Outputs:
%   star_meas - Cell array (size N) or 3D array (3 x M x N)
%               Here we assume all stars are visible at all times for simplicity.
%               Output is 3 x M x N matrix.

    [~, N] = size(q_true);
    [~, M] = size(v_ref);
    
    star_meas = zeros(3, M, N);
    
    sigma_st = params.st_noise_std;
    
    for k = 1:N
        q_k = q_true(:, k);
        A_k = utils_AttitudeMatrix(q_k);
        
        for i = 1:M
            curr_ref = v_ref(:, i);
            
            % True Body Vector
            v_body_true = A_k * curr_ref;
            
            % Add Noise
            % Simple approach: Add additive Gaussian noise and re-normalize
            noise = sigma_st * randn(3, 1);
            v_meas_raw = v_body_true + noise;
            
            % Normalize
            v_meas = v_meas_raw / norm(v_meas_raw);
            
            star_meas(:, i, k) = v_meas;
        end
    end

end
