% MAIN_SIMULATION  PPV-aided SWGLT fault detection for INS/ISIS redundant system.
%   Generates flight trajectory, simulates sensor measurements for 2 INS + 1 ISIS,
%   injects faults under three conditions, runs three FDI algorithms (GLT, WGLT,
%   PPV-aided SWGLT), evaluates performance, and plots results.
%
% Requires: PSINS toolbox initialized (run psinsinit.m from PSINS root first).
%
% See also  config, inject_fault, fdi_glt, fdi_wglt, fdi_swglt,
%           evaluate_performance, plot_results.

clear; close all; clc;

%% Initialize PSINS toolbox
glvs;
rng(42);

%% Load configuration
cfg = config();
ts = cfg.ts;

%% ========================================================================
%  Step 1: Generate flight trajectory using NLG approach
% =========================================================================
fprintf('=== Step 1: Generating flight trajectory ===\n');

trjPath = fullfile(pwd, 'data', 'trj_NLG_approach.mat');
if exist(trjPath, 'file')
    trj = trjfile(trjPath);
    fprintf('Loaded existing trajectory from %s\n', trjPath);
else
    trj = generate_trajectory(ts);
    trjfile(trjPath, trj);
    fprintf('Generated and saved trajectory to %s\n', trjPath);
end

% Extract reference IMU data
imu_ref = trj.imu;          % Nx7: [dtheta_xyz, dv_xyz, time]
N = size(imu_ref, 1);
t = imu_ref(:, 7);          % time vector (s)
T_total = t(end) - t(1);

fprintf('Trajectory: %.1f s, %d samples, ts=%.3f s\n', T_total, N, ts);

%% ========================================================================
%  Step 2: Generate redundant sensor measurements (2 INS + 1 ISIS)
% =========================================================================
fprintf('\n=== Step 2: Generating sensor measurements ===\n');

% Convert angular/velocity increments to rates/specific forces
gyro_ref = imu_ref(:, 1:3) / ts;    % angular rate (rad/s)
acc_ref  = imu_ref(:, 4:6) / ts;    % specific force (m/s^2)

% Generate noise for each subsystem
% INS1 (sensors 1-3)
noise_ins1_gyro = cfg.sigma_ins_gyro  * randn(N, 3);
noise_ins1_acc  = cfg.sigma_ins_acc   * randn(N, 3);
% INS2 (sensors 4-6)
noise_ins2_gyro = cfg.sigma_ins_gyro  * randn(N, 3);
noise_ins2_acc  = cfg.sigma_ins_acc   * randn(N, 3);
% ISIS (sensors 7-9)
noise_isis_gyro = cfg.sigma_isis_gyro * randn(N, 3);
noise_isis_acc  = cfg.sigma_isis_acc  * randn(N, 3);

% Assemble 9-sensor measurement vectors: Z = H * x_true + noise
Z_gyro = [gyro_ref + noise_ins1_gyro, ...
          gyro_ref + noise_ins2_gyro, ...
          gyro_ref + noise_isis_gyro];       % N x 9

Z_acc  = [acc_ref  + noise_ins1_acc, ...
          acc_ref  + noise_ins2_acc, ...
          acc_ref  + noise_isis_acc];        % N x 9

fprintf('Generated %d samples x 9 sensors (gyro + acc)\n', N);

%% ========================================================================
%  Step 3: Run FDI for each fault condition
% =========================================================================
conditions = {'Condition 1: INS1 Y-axis Gyro Hard Fault', ...
              'Condition 2: INS2 X-axis Accel Hard Fault', ...
              'Condition 3: INS1 Z-axis Gyro Soft Fault'};

for cond = 1:3
    fprintf('\n============================================================\n');
    fprintf('=== %s ===\n', conditions{cond});
    fprintf('============================================================\n');

    %% Inject faults
    switch cond
        case 1
            Z_work     = inject_fault(Z_gyro, t, cfg.fault1, 'hard');
            sigma_vec  = cfg.sigma_gyro;
            sigma_uni  = cfg.sigma_ins_gyro;   % uniform sigma for GLT
            fault_cfg  = cfg.fault1;
        case 2
            Z_work     = inject_fault(Z_acc, t, cfg.fault2, 'hard');
            sigma_vec  = cfg.sigma_acc;
            sigma_uni  = cfg.sigma_ins_acc;
            fault_cfg  = cfg.fault2;
        case 3
            Z_work     = inject_fault(Z_gyro, t, cfg.fault3, 'soft');
            sigma_vec  = cfg.sigma_gyro;
            sigma_uni  = cfg.sigma_ins_gyro;
            fault_cfg  = cfg.fault3;
            if isfield(fault_cfg, 't_start_ref') && ~isempty(fault_cfg.t_start_ref)
                t_ref = fault_cfg.t_start_ref;
            else
                t_ref = fault_cfg.interval(1);
            end
            fprintf('Soft fault config: interval=[%.1f, %.1f] s, rate=%.4e, t_ref=%.1f s\n', ...
                fault_cfg.interval(1), fault_cfg.interval(2), fault_cfg.rate, t_ref);
    end

    % Build fault time mask
    fault_mask = false(N, 1);
    if cond <= 2
        for fi = 1:size(fault_cfg.intervals, 1)
            fault_mask = fault_mask | ...
                (t >= fault_cfg.intervals(fi,1) & t <= fault_cfg.intervals(fi,2));
        end
    else
        fault_mask = (t >= fault_cfg.interval(1) & t <= fault_cfg.interval(2));
    end

    %% Run three algorithms
    fprintf('\n--- Running GLT ---\n');
    [FD_glt,  FI_glt,  V_glt]  = fdi_glt(Z_work, sigma_uni, cfg);

    fprintf('--- Running WGLT ---\n');
    [FD_wglt, FI_wglt, V_wglt] = fdi_wglt(Z_work, sigma_vec, cfg);

    fprintf('--- Running SWGLT ---\n');
    [FD_swglt, FI_swglt, T_adaptive] = fdi_swglt(Z_work, sigma_vec, cfg);

    %% Evaluate performance
    fprintf('\n--- Performance Evaluation ---\n');
    is_soft = (cond == 3);

    stats_glt   = evaluate_performance(FD_glt,   cfg.T_D,     fault_mask, 'GLT',   is_soft, t);
    stats_wglt  = evaluate_performance(FD_wglt,  cfg.T_D,     fault_mask, 'WGLT',  is_soft, t);
    stats_swglt = evaluate_performance(FD_swglt, T_adaptive,  fault_mask, 'SWGLT', is_soft, t);

    % Summary table
    fprintf('\n  Summary:\n');
    fprintf('  %-8s  FDR=%6.2f%%  FAR=%6.2f%%  Accuracy=%6.2f%%', ...
        'GLT', stats_glt.FDR*100, stats_glt.FAR*100, stats_glt.Accuracy*100);
    if is_soft, fprintf('  Delay=%.1fs', stats_glt.delay); end; fprintf('\n');

    fprintf('  %-8s  FDR=%6.2f%%  FAR=%6.2f%%  Accuracy=%6.2f%%', ...
        'WGLT', stats_wglt.FDR*100, stats_wglt.FAR*100, stats_wglt.Accuracy*100);
    if is_soft, fprintf('  Delay=%.1fs', stats_wglt.delay); end; fprintf('\n');

    fprintf('  %-8s  FDR=%6.2f%%  FAR=%6.2f%%  Accuracy=%6.2f%%', ...
        'SWGLT', stats_swglt.FDR*100, stats_swglt.FAR*100, stats_swglt.Accuracy*100);
    if is_soft, fprintf('  Delay=%.1fs', stats_swglt.delay); end; fprintf('\n');

    %% Plot results
    if cond <= 2
        fault_intervals = fault_cfg.intervals;
    else
        fault_intervals = fault_cfg.interval;
    end

    plot_results(t, FD_glt, FD_wglt, FD_swglt, ...
                 cfg.T_D, T_adaptive, ...
                 FI_glt, FI_wglt, FI_swglt, ...
                 fault_mask, fault_intervals, ...
                 cond, conditions{cond}, cfg);
end

fprintf('\n=== Simulation Complete ===\n');


%% ========================================================================
%  Local function: Generate trajectory
% =========================================================================
function trj = generate_trajectory(ts)
    global glv

    legSpeeds = [210, 200, 190, 170, 150, 140] * glv.kn;
    accelTime = 60;  decelTime = 40;
    turnRate = 3;    turnRollTime = 6;
    pitchRate = 0.5;

    % Waypoint definition
    W = struct('name',{},'role',{},'lat',{},'lon',{},'alt',{},'alt_src',{},'E',{},'N',{});
    W(end+1) = makeWP('NLG',   'IAF',    22.5317, 113.5617, 4900, 'MSL');
    W(end+1) = makeWP('SZ928', 'Fly-by', 22.4867, 113.6517, 2600, 'MSL');
    W(end+1) = makeWP('SZA34', 'Fly-by', 22.4417, 113.7400, 2300, 'MSL');
    W(end+1) = makeWP('SZ924', 'IF',     22.4683, 113.7967, 1300, 'MSL');
    W(end+1) = makeWP('SZ923', 'Fly-by', 22.4783, 113.8183, 1100, 'MSL');
    W(end+1) = makeWP('SZ921', 'Fly-by', 22.4867, 113.8350, 1100, 'MSL');
    W(end+1) = makeWP('SZ920', 'RF_End', 22.5333, 113.8517, 1100, 'MSL');
    W = fillWaypointAltitude(W);

    legs = waypointLegs(W, legSpeeds);
    avp0 = [[0; 0; legs(1).course]; [0; 0; 0]; [W(1).lat; W(1).lon; W(1).alt]];

    seg = trjsegment([], 'init', 0);
    speedAccelRate = legSpeeds(1) / accelTime;
    speedDecelRate = legSpeeds(1) / decelTime;
    currentSpeed = 0;
    targetSpeed = legSpeeds(1);
    seg = adjustSpeed(seg, targetSpeed, currentSpeed, speedAccelRate, speedDecelRate);
    currentSpeed = targetSpeed;

    currentPitch = 0;
    currentYaw = legs(1).course;
    if abs(legs(1).pitch) > 0.1*glv.deg
        pitchTime = abs(legs(1).pitch) / glv.deg / pitchRate;
        if legs(1).pitch > 0
            seg = trjsegment(seg, 'headup', pitchTime, pitchRate);
        else
            seg = trjsegment(seg, 'headdown', pitchTime, pitchRate);
        end
        currentPitch = legs(1).pitch;
    end
    seg = trjsegment(seg, 'uniform', legs(1).time);

    for k = 2:numel(legs)
        deltaYaw = wrapToPi_(legs(k).course - currentYaw);
        deltaYawDeg = deltaYaw / glv.deg;
        if abs(deltaYawDeg) > 0.5
            turnTime = abs(deltaYawDeg) / turnRate;
            if deltaYawDeg > 0
                seg = trjsegment(seg, 'coturnleft', turnTime, turnRate, [], turnRollTime);
            else
                seg = trjsegment(seg, 'coturnright', turnTime, turnRate, [], turnRollTime);
            end
            currentYaw = currentYaw + deltaYaw;
        end
        deltaPitch = legs(k).pitch - currentPitch;
        deltaPitchDeg = deltaPitch / glv.deg;
        if abs(deltaPitchDeg) > 0.1
            pitchTime = abs(deltaPitchDeg) / pitchRate;
            if deltaPitchDeg > 0
                seg = trjsegment(seg, 'headup', pitchTime, pitchRate);
            else
                seg = trjsegment(seg, 'headdown', pitchTime, pitchRate);
            end
            currentPitch = currentPitch + deltaPitch;
        end
        targetSpeed = legs(k).speed;
        seg = adjustSpeed(seg, targetSpeed, currentSpeed, speedAccelRate, speedDecelRate);
        currentSpeed = targetSpeed;
        seg = trjsegment(seg, 'uniform', legs(k).time);
    end

    if abs(currentPitch) > 0.1*glv.deg
        pitchTime = abs(currentPitch) / glv.deg / pitchRate;
        if currentPitch > 0
            seg = trjsegment(seg, 'headdown', pitchTime, pitchRate);
        else
            seg = trjsegment(seg, 'headup', pitchTime, pitchRate);
        end
    end
    seg = adjustSpeed(seg, 0, currentSpeed, speedAccelRate, speedDecelRate);
    seg = trjsegment(seg, 'uniform', 20);

    % Generate at 0.01s then downsample to cfg.ts
    trj_raw = trjsimu(avp0, seg.wat, 0.01, 1);

    if abs(ts - 0.01) > 1e-6
        ratio = round(ts / 0.01);
        nSamp = floor(size(trj_raw.imu, 1) / ratio);
        trj.imu = zeros(nSamp, 7);
        for idx = 1:nSamp
            rows = ((idx-1)*ratio+1) : (idx*ratio);
            trj.imu(idx, 1:3) = sum(trj_raw.imu(rows, 1:3), 1);
            trj.imu(idx, 4:6) = sum(trj_raw.imu(rows, 4:6), 1);
            trj.imu(idx, 7)   = trj_raw.imu(rows(end), 7);
        end
        trj.avp = trj_raw.avp(ratio:ratio:end, :);
        if size(trj.avp, 1) > nSamp
            trj.avp = trj.avp(1:nSamp, :);
        end
    else
        trj.imu = trj_raw.imu;
        trj.avp = trj_raw.avp;
    end
    trj.avp0 = trj_raw.avp0;
    trj.ts = ts;
end

function wp = makeWP(name, role, lat_deg, lon_deg, alt, alt_src)
    global glv
    wp.name = name;  wp.role = role;
    wp.lat = lat_deg * glv.deg;  wp.lon = lon_deg * glv.deg;
    wp.alt = alt;  wp.alt_src = alt_src;
    wp.E = NaN;  wp.N = NaN;
end

function W = fillWaypointAltitude(W)
    alts = [W.alt]';
    idx = find(~isnan(alts));
    if numel(idx) >= 2
        alts = interp1(idx, alts(idx), 1:numel(W), 'linear', 'extrap')';
    else
        alts(isnan(alts)) = alts(idx);
    end
    for k = 1:numel(W), W(k).alt = alts(k); end
end

function legs = waypointLegs(W, speed)
    global glv
    n = numel(W);
    legs = repmat(struct('label','','course',0,'pitch',0,...
        'horizontal',0,'vertical',0,'track',0,'time',0,'speed',0), n-1, 1);
    speeds = speed(:);
    for k = 1:n-1
        pos1 = [W(k).lat, W(k).lon, W(k).alt];
        pos2 = [W(k+1).lat, W(k+1).lon, W(k+1).alt];
        dxyz = pos2dxyz([pos1; pos2], pos1');
        enu = dxyz(end, 1:3);
        horizontal = hypot(enu(1), enu(2));
        vertical = enu(3);
        legs(k).label = sprintf('%s -> %s', W(k).name, W(k+1).name);
        legs(k).course = wrapToPi_(atan2(enu(1), enu(2)));
        legs(k).pitch = atan2(vertical, horizontal);
        legs(k).horizontal = horizontal;
        legs(k).vertical = vertical;
        legs(k).track = sqrt(horizontal^2 + vertical^2);
        legs(k).speed = speeds(k);
        legs(k).time = horizontal / speeds(k);
    end
end

function ang = wrapToPi_(ang)
    ang = mod(ang + pi, 2*pi) - pi;
end

function seg = adjustSpeed(seg, targetSpeed, currentSpeed, accelRate, decelRate)
    deltaV = targetSpeed - currentSpeed;
    if deltaV > 1e-3
        seg = trjsegment(seg, 'accelerate', deltaV / accelRate, [], accelRate);
    elseif deltaV < -1e-3
        seg = trjsegment(seg, 'deaccelerate', -deltaV / decelRate, [], decelRate);
    end
end
