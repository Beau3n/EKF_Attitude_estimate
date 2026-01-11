function [gyro_meas, b_true] = Model_Gyro_Sensor(w_true, dt, params)
% MODEL_GYRO_SENSOR Simulates Gyroscope Measurements with Bias and Noise
%
% Model: u_meas = w_true + b + eta_v
%        db/dt = eta_u (Bias Random Walk)
%
% Inputs:
%   w_true - True angular velocity (3 x N)
%   dt     - Sampling time
%   params - Struct with fields:
%            params.gyro_noise_std  (ARW: Angle Random Walk)
%            params.gyro_bias_std   (RRW: Rate Random Walk / Bias Instability)
%            params.b_init          (Initial bias)
%
% Outputs:
%   gyro_meas - Measured angular velocity (3 x N)
%   b_true    - True bias history (3 x N)

    [~, N] = size(w_true);
    
    gyro_meas = zeros(3, N);
    b_true = zeros(3, N);
    
    % Noise parameters
    sigma_v = params.gyro_noise_std; % rad/s
    sigma_u = params.gyro_bias_std;  % rad/s/sqrt(s) -> Discrete update scaling needed
    
    % Bias Random Walk scaling: sigma_bd = sigma_u * sqrt(dt)
    sigma_bd = sigma_u * sqrt(dt);
    
    % Initial Bias
    b_curr = params.b_init;
    
    rng(42); % Seed for reproducibility
    
    for k = 1:N
        % Store true bias
        b_true(:, k) = b_curr;
        
        % Measurement Noise (White Noise)
        eta_v = sigma_v * randn(3, 1);
        
        % Measurement
        gyro_meas(:, k) = w_true(:, k) + b_curr + eta_v;
        
        % Bias Random Walk Update
        eta_u = sigma_bd * randn(3, 1);
        b_curr = b_curr + eta_u;
    end

end
