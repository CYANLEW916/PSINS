% Sliding-window whitened GLRT FDI demo for a heterogeneous 2 INS + 1 ISIS architecture.
% Inputs:  data/trj_NLG_approach.mat (trj.imu as reference IMU increment sequence).
% Outputs: detection/isolation plots and adaptive threshold diagnostics.

glvs;
close all;

%% Step 1: Load trajectory IMU reference data
trjData = load('data/trj_NLG_approach.mat');
trj = trjData.trj;
imu_ref = trj.imu;               % Nx7: [dtheta_xyz, dv_xyz, time]
N = size(imu_ref, 1);
t = imu_ref(:, 7);
ts = mean(diff(t));
T_total = t(end) - t(1);

fprintf('Trajectory: %.1f s, %d samples, ts=%.3f s\n', T_total, N, ts);

gyro_ref = imu_ref(:, 1:3) / ts; % angular rate (rad/s)
acc_ref = imu_ref(:, 4:6) / ts;  % specific force (m/s^2)

%% Step 2: Build heterogeneous redundant model (2 INS + 1 ISIS)
n = 9; m = 3;
H = [eye(3); eye(3); eye(3)];

sigma_ins_gyro = 0.01 * pi / (180 * 3600);
sigma_isis_gyro = 0.1 * pi / (180 * 3600);
Rgyro = diag([repmat(sigma_ins_gyro^2, 1, 6), repmat(sigma_isis_gyro^2, 1, 3)]);

j = 25;
K = 0.01;
PdStar = 0.90;
modelGyro = fdi_precompute_model(Rgyro, H, K, PdStar, j);

fprintf('Adaptive threshold (gyro): %.4f, fixed threshold (Pfa=0.01): %.4f\n', ...
    modelGyro.Tadapt, modelGyro.Td);
fprintf('Predicted adaptive Pfa: %.6f\n', modelGyro.PfaAdapt);

%% Step 3: Generate redundant gyro measurements and inject a fault
faultCfg.sensorIdx = 8;
faultCfg.startIdx = round(0.55 * N);
faultCfg.type = 'step';
faultCfg.amplitude = 3.0e-5;  % rad/s

[Zgyro, gyroMeta] = fdi_generate_redundant_measurements(
    gyro_ref, sigma_ins_gyro, sigma_isis_gyro, faultCfg);

%% Step 4: Run sliding-window WGLT detection/isolation
out = fdi_run_sliding_wglt(Zgyro, modelGyro);
firstAlarm = find(out.alarm, 1, 'first');

if isempty(firstAlarm)
    fprintf('No fault alarm was triggered.\n');
else
    fprintf('First alarm at t=%.2f s, isolated sensor=%d, fhat=%.3e rad/s\n', ...
        t(firstAlarm), out.kstar(firstAlarm), out.fhat(firstAlarm));
end

%% Step 5: Visualization
figure('Name', 'WGLT FDI demo: gyro');
subplot(3, 1, 1);
plot(t, out.FDj, 'b', 'LineWidth', 1.0); hold on;
yline(modelGyro.Tadapt, 'r--', 'LineWidth', 1.0);
grid on; ylabel('FD_j');
title('Sliding-window energy statistic and adaptive threshold');
legend('FD_j', 'T_{adapt}', 'Location', 'northwest');

subplot(3, 1, 2);
plot(t, out.kstar, 'k', 'LineWidth', 1.0); hold on;
yline(faultCfg.sensorIdx, 'r--', 'LineWidth', 1.0);
grid on; ylabel('isolated idx');
title('Isolation result');
legend('k^*', 'true fault idx', 'Location', 'northwest');

subplot(3, 1, 3);
plot(t, out.fhat, 'b', 'LineWidth', 1.0); hold on;
yline(faultCfg.amplitude, 'r--', 'LineWidth', 1.0);
grid on; ylabel('fhat (rad/s)'); xlabel('time (s)');
title('Fault magnitude estimate');
legend('fhat', 'true fault', 'Location', 'northwest');

%% Optional: same pipeline can be repeated for acc_ref with accel noise levels.
% sigma_ins_acc = 1e-4 * 9.80665;
% sigma_isis_acc = 5e-3 * 9.80665;
% Racc = diag([repmat(sigma_ins_acc^2, 1, 6), repmat(sigma_isis_acc^2, 1, 3)]);
% modelAcc = fdi_precompute_model(Racc, H, K, PdStar, j);
% [Zacc, ~] = fdi_generate_redundant_measurements(acc_ref, sigma_ins_acc, sigma_isis_acc);
% outAcc = fdi_run_sliding_wglt(Zacc, modelAcc);
