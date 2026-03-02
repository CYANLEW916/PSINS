function summaryTbl = test_FDI_WGLT_fault_campaign()
%TEST_FDI_WGLT_FAULT_CAMPAIGN Evaluate WGLT under configured fault campaign.
% Fault campaign evaluation for heterogeneous 2 INS + 1 ISIS WGLT FDI.
% Inputs:
%   data/trj_NLG_approach.mat (trj.imu).
% Outputs:
%   Per-case visualization and summary metrics: detection rate, isolation accuracy,
%   and first-alarm delay for each fault window.

glvs;
close all;
rng(2026);

%% Step 1: Load trajectory and extract reference signals
trjData = load('data/trj_NLG_approach.mat');
trj = trjData.trj;
imu_ref = trj.imu;                    % [dtheta_xyz, dv_xyz, t]
N = size(imu_ref, 1);
t = imu_ref(:, 7);
ts = mean(diff(t));

gyro_ref = imu_ref(:, 1:3) / ts;      % rad/s
acc_ref = imu_ref(:, 4:6) / ts;       % m/s^2
fprintf('Trajectory loaded: N=%d, ts=%.3f s, T=%.1f s\n', N, ts, t(end)-t(1));

%% Step 2: Model setup
H = [eye(3); eye(3); eye(3)];
j = 25;
K = 0.01;
PdStar = 0.90;

sigma_ins_gyro = 0.01 * pi / (180 * 3600);
sigma_isis_gyro = 0.1 * pi / (180 * 3600);
Rgyro = diag([repmat(sigma_ins_gyro^2, 1, 6), repmat(sigma_isis_gyro^2, 1, 3)]);
modelGyro = fdi_precompute_model(Rgyro, H, K, PdStar, j);

sigma_ins_acc = 1e-4 * 9.80665;
sigma_isis_acc = 5e-3 * 9.80665;
Racc = diag([repmat(sigma_ins_acc^2, 1, 6), repmat(sigma_isis_acc^2, 1, 3)]);
modelAcc = fdi_precompute_model(Racc, H, K, PdStar, j);

%% Step 3: Define fault conditions from paper table
% Channel mapping for [INS1(x,y,z), INS2(x,y,z), ISIS(x,y,z)]
cases = {};

c1.name = 'Case 1: IRS1 Y-gyro step faults';
c1.sensorType = 'gyro';
c1.events = [...
    struct('sensorIdx', 2, 'startT', 120, 'endT', 135, 'type', 'step', 'amp', 0.5*pi/(180*3600)); ...
    struct('sensorIdx', 2, 'startT', 440, 'endT', 455, 'type', 'step', 'amp', 1.0*pi/(180*3600)); ...
    struct('sensorIdx', 2, 'startT', 595, 'endT', 610, 'type', 'step', 'amp', 2.0*pi/(180*3600)) ...
    ];
cases{end+1} = c1;

c2.name = 'Case 2: IRS2 X-accel step faults';
c2.sensorType = 'acc';
c2.events = [...
    struct('sensorIdx', 4, 'startT', 100, 'endT', 115, 'type', 'step', 'amp', 0.005*9.80665); ...
    struct('sensorIdx', 4, 'startT', 420, 'endT', 435, 'type', 'step', 'amp', 0.01*9.80665); ...
    struct('sensorIdx', 4, 'startT', 575, 'endT', 590, 'type', 'step', 'amp', 0.02*9.80665) ...
    ];
cases{end+1} = c2;

c3.name = 'Case 3: IRS1 Z-gyro ramp fault';
c3.sensorType = 'gyro';
c3.events = struct('sensorIdx', 3, 'startT', 260, 'endT', 310, ...
    'type', 'ramp', 'amp', 0.02*pi/(180*3600));
cases{end+1} = c3;

%% Step 4: Run campaign and evaluate
summary = cell(numel(cases), 6);

for ic = 1:numel(cases)
    ci = cases{ic};
    evtIdx = time_events_to_index(ci.events, t);

    faultCfg.events = evtIdx;
    if strcmpi(ci.sensorType, 'gyro')
        [Z, meta] = fdi_generate_redundant_measurements(
            gyro_ref, sigma_ins_gyro, sigma_isis_gyro, faultCfg);
        out = fdi_run_sliding_wglt(Z, modelGyro);
        threshold = modelGyro.Tadapt;
        unitStr = 'rad/s';
    else
        [Z, meta] = fdi_generate_redundant_measurements(
            acc_ref, sigma_ins_acc, sigma_isis_acc, faultCfg);
        out = fdi_run_sliding_wglt(Z, modelAcc);
        threshold = modelAcc.Tadapt;
        unitStr = 'm/s^2';
    end

    metrics = evaluate_case(out, meta.fault, evtIdx, t);
    summary{ic, 1} = ci.name;
    summary{ic, 2} = metrics.totalWindows;
    summary{ic, 3} = metrics.detectedWindows;
    summary{ic, 4} = metrics.detectionRate;
    summary{ic, 5} = metrics.isolationAccuracy;
    summary{ic, 6} = metrics.meanDelay;

    plot_case(ci.name, t, out, threshold, meta.fault, evtIdx, unitStr);

    fprintf('\n%s\n', ci.name);
    fprintf('  windows: %d, detected: %d, detection rate: %.2f%%\n', ...
        metrics.totalWindows, metrics.detectedWindows, 100*metrics.detectionRate);
    fprintf('  isolation accuracy: %.2f%%, mean delay: %.2f s\n', ...
        100*metrics.isolationAccuracy, metrics.meanDelay);
end

summaryTbl = cell2table(summary, 'VariableNames', ...
    {'CaseName', 'Windows', 'Detected', 'DetectionRate', 'IsolationAcc', 'MeanDelaySec'});
disp('=== WGLT fault campaign summary ===');
disp(summaryTbl);

end

function evtIdx = time_events_to_index(events, t)
%TIME_EVENTS_TO_INDEX Convert event windows from seconds to sample indices.

    nEvt = numel(events);
    evtIdx = repmat(struct('sensorIdx', 0, 'startIdx', 1, 'endIdx', 1, ...
        'type', 'step', 'amplitude', 0), nEvt, 1);

    for k = 1:nEvt
        evtIdx(k).sensorIdx = events(k).sensorIdx;
        k0 = find(t >= events(k).startT, 1, 'first');
        k1 = find(t <= events(k).endT, 1, 'last');
        if isempty(k0), k0 = 1; end
        if isempty(k1), k1 = numel(t); end
        evtIdx(k).startIdx = k0;
        evtIdx(k).endIdx = k1;
        evtIdx(k).type = events(k).type;
        evtIdx(k).amplitude = events(k).amp;
    end
end

function metrics = evaluate_case(out, fault, events, t)
%EVALUATE_CASE Compute detection/isolation metrics for each fault window.

    totalWindows = numel(events);
    detectedWindows = 0;
    correctIsolation = 0;
    delaySum = 0;

    for k = 1:totalWindows
        k0 = events(k).startIdx;
        k1 = events(k).endIdx;
        idxWin = k0:k1;
        idxAlarm = idxWin(out.alarm(idxWin));

        if ~isempty(idxAlarm)
            detectedWindows = detectedWindows + 1;
            delaySum = delaySum + (t(idxAlarm(1)) - t(k0));

            if any(out.kstar(idxAlarm) == events(k).sensorIdx)
                correctIsolation = correctIsolation + 1;
            end
        end
    end

    if detectedWindows > 0
        meanDelay = delaySum / detectedWindows;
        isoAcc = correctIsolation / detectedWindows;
    else
        meanDelay = NaN;
        isoAcc = 0;
    end

    metrics.totalWindows = totalWindows;
    metrics.detectedWindows = detectedWindows;
    metrics.detectionRate = detectedWindows / totalWindows;
    metrics.isolationAccuracy = isoAcc;
    metrics.meanDelay = meanDelay;
    metrics.faultOccupancy = nnz(any(fault ~= 0, 2)) / size(fault, 1);
end

function plot_case(nameStr, t, out, threshold, fault, events, unitStr)
%PLOT_CASE Visualize FD statistic, isolation channel, and true/estimated fault traces.

    figure('Name', nameStr);

    subplot(3, 1, 1);
    plot(t, out.FDj, 'b', 'LineWidth', 1.0); hold on;
    yline(threshold, 'r--', 'LineWidth', 1.0);
    for k = 1:numel(events)
        xline(t(events(k).startIdx), 'k:');
        xline(t(events(k).endIdx), 'k:');
    end
    grid on;
    ylabel('FD_j');
    title([nameStr, ' - Detection']);
    legend('FD_j', 'T_{adapt}', 'Location', 'northwest');

    subplot(3, 1, 2);
    plot(t, out.kstar, 'k', 'LineWidth', 1.0); hold on;
    for k = 1:numel(events)
        yline(events(k).sensorIdx, 'r--');
    end
    grid on;
    ylabel('isolated idx');
    title('Isolation index');

    subplot(3, 1, 3);
    trueFault = sum(fault, 2);
    plot(t, trueFault, 'r', 'LineWidth', 1.0); hold on;
    plot(t, out.fhat, 'b', 'LineWidth', 1.0);
    grid on;
    ylabel(['fault (', unitStr, ')']);
    xlabel('time (s)');
    title('True vs estimated fault magnitude');
    legend('true fault (aggregated)', 'estimated f_s', 'Location', 'northwest');
end
