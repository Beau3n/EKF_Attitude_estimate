% Test_Simulation_Setup.m
% A temporary script to verify the simulation environment and utils.

clc; clear; close all;
addpath('Utils');
addpath('Simulation');

% 1. Config
dt = 0.1;
T_total = 100;
params.gyro_noise_std = 1e-4;
params.gyro_bias_std = 1e-5;
params.b_init = [0.01; -0.01; 0.005];
params.st_noise_std = 1e-3;

% 2. Kinematics
fprintf('Running Kinematics Simulation...\n');
[time, q_true, w_true] = sim_Kinematics(T_total, dt);

% Check Quaternion Norm
norms = sqrt(sum(q_true.^2, 1));
fprintf('Max Latent Quaternion Norm Error: %e\n', max(abs(norms - 1)));

% 3. Sensor Models
fprintf('Generating Sensor Data...\n');
[gyro_meas, b_true] = Model_Gyro_Sensor(w_true, dt, params);

% Define 2 Reference Stars
v_ref = [1, 0; 0, 1; 0, 0]; % Star 1 x-axis, Star 2 y-axis
star_meas = Model_StarTracker(q_true, v_ref, params);

% 4. Visualization
figure;
subplot(3,1,1);
plot(time, q_true); title('True Quaternion'); legend('q1','q2','q3','q4');
grid on;

subplot(3,1,2);
plot(time, w_true, '--'); hold on;
plot(time, gyro_meas, '-'); 
title('Angular Velocity: True vs Measured');
legend('wx','wy','wz');
grid on;

subplot(3,1,3);
plot(time, b_true); title('True Bias');
grid on;

fprintf('Simulation Test Complete.\n');
