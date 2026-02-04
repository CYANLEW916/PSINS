% Fault injection and detection evaluation for IRS1/2 and ISIS sensors.
% Requires 'trj10ms_sensor_data.mat' from test_SINS_IRS_ISIS.m before running.
% See also  test_SINS_IRS_ISIS, test_SINS_trj.

glvs

if exist('trj10ms_sensor_data.mat', 'file') ~= 2
    test_SINS_IRS_ISIS;
end

load('trj10ms_sensor_data.mat', 'sensorData', 'specs', 'trj');

faultIdx = find(strcmp({sensorData.name}, 'ISIS'), 1, 'first');
if isempty(faultIdx)
    error('ISIS sensor data not found.');
end

avpNom = sensorData(faultIdx).avp;
t = avpNom(:, end);
trjAvp = trj.avp(1:size(avpNom, 1), :);

faultCfg = buildFaultConfig(t);
detCfg = makeDetectionConfig();

[avpFault, faultMask] = injectFaults(avpNom, t, faultCfg);
[detMask, resFault] = detectFaults(avpFault, trjAvp, t, detCfg);
[detMaskNom, resNom] = detectFaults(avpNom, trjAvp, t, detCfg);

metricsFault = calcFdiMetrics(faultMask, detMask, t);
metricsNom = calcFdiMetrics(false(size(t)), detMaskNom, t);

avpMit = avpFault;
avpMit(detMask, :) = avpNom(detMask, :);

rmseNom = computeRmse(avpNom, trjAvp);
rmseFault = computeRmse(avpFault, trjAvp);
rmseMit = computeRmse(avpMit, trjAvp);

fprintf(['\nFault detection metrics (faulty): FDR=%.3f, FAR=%.3f, ', ...
    'delay=%.2fs\n'], metricsFault.fdr, metricsFault.far, metricsFault.delay);
fprintf(['Fault detection metrics (nominal): FDR=%.3f, FAR=%.3f, ', ...
    'delay=%.2fs\n'], metricsNom.fdr, metricsNom.far, metricsNom.delay);
fprintf(['RMSE att/vel/pos (nominal): [%.3f %.3f %.3f] ', ...
    '[%.3f %.3f %.3f] [%.3f %.3f %.3f]\n'], rmseNom.att, rmseNom.vel, ...
    rmseNom.pos);
fprintf(['RMSE att/vel/pos (faulty): [%.3f %.3f %.3f] ', ...
    '[%.3f %.3f %.3f] [%.3f %.3f %.3f]\n'], rmseFault.att, rmseFault.vel, ...
    rmseFault.pos);
fprintf(['RMSE att/vel/pos (mitigated): [%.3f %.3f %.3f] ', ...
    '[%.3f %.3f %.3f] [%.3f %.3f %.3f]\n'], rmseMit.att, rmseMit.vel, ...
    rmseMit.pos);

save('trj10ms_sensor_data_faults.mat', 'specs', 'faultCfg', 'detCfg', ...
    'metricsFault', 'metricsNom', 'rmseNom', 'rmseFault', 'rmseMit', ...
    'avpFault', 'avpMit', 't');

plotFdiResults(t, resFault, detMask, faultMask, detCfg, 'Faulty case');
plotFdiResults(t, resNom, detMaskNom, false(size(t)), detCfg, ...
    'Nominal case');
plotNavErrors(t, avpNom, avpFault, avpMit, trjAvp);

%% helper functions
function faultCfg = buildFaultConfig(t)
%BUILDFAULTCONFIG Configure detectable faults for the ISIS sensor.
%   faultCfg = BUILDFAULTCONFIG(t) returns a struct array with fault windows
%   and magnitudes based on the time vector t.
    tEnd = t(end);
    faultCfg = struct( ...
        'name', {'att_step', 'vel_step', 'pos_ramp'}, ...
        'tStart', {0.2 * tEnd, 0.5 * tEnd, 0.7 * tEnd}, ...
        'tEnd', {0.3 * tEnd, 0.55 * tEnd, 0.85 * tEnd}, ...
        'attDeg', {[2 0 0], [], []}, ...
        'vel', {[], [0.5 0 0], []}, ...
        'posRampM', {[], [], [50 0 0]});
end

function detCfg = makeDetectionConfig()
%MAKEDETECTIONCONFIG Thresholds for fault detection.
%   detCfg = MAKEDETECTIONCONFIG() returns thresholds and persistence.
    detCfg = struct('attDeg', 0.6, 'vel', 0.3, 'posM', 15, ...
        'minPersist', 5);
end

function [avpFault, faultMask] = injectFaults(avp, t, faultCfg)
%INJECTFAULTS Inject additive faults into AVP outputs.
%   [avpFault, faultMask] = INJECTFAULTS(avp, t, faultCfg) injects faults
%   defined in faultCfg into avp and returns the fault mask.
    global glv
    avpFault = avp;
    faultMask = false(size(t));
    for k = 1:numel(faultCfg)
        idx = t >= faultCfg(k).tStart & t <= faultCfg(k).tEnd;
        faultMask = faultMask | idx;
        if ~isempty(faultCfg(k).attDeg)
            attAdd = repmat(faultCfg(k).attDeg, sum(idx), 1) * glv.deg;
            avpFault(idx, 1:3) = avpFault(idx, 1:3) + attAdd;
        end
        if ~isempty(faultCfg(k).vel)
            velAdd = repmat(faultCfg(k).vel, sum(idx), 1);
            avpFault(idx, 4:6) = avpFault(idx, 4:6) + velAdd;
        end
        if ~isempty(faultCfg(k).posRampM)
            ramp = (t(idx) - faultCfg(k).tStart) / ...
                (faultCfg(k).tEnd - faultCfg(k).tStart);
            ramp = ramp(:);
            dpos = bsxfun(@times, ramp, faultCfg(k).posRampM);
            dpos(:, 1:2) = dpos(:, 1:2) / glv.Re;
            avpFault(idx, 7:9) = avpFault(idx, 7:9) + dpos;
        end
    end
end

function [detMask, res] = detectFaults(avp, avpRef, t, detCfg)
%DETECTFAULTS Detect faults based on residual thresholds.
%   [detMask, res] = DETECTFAULTS(avp, avpRef, t, detCfg) returns the
%   detection mask and residual diagnostics.
    global glv
    attErrDeg = (avp(:, 1:3) - avpRef(:, 1:3)) / glv.deg;
    velErr = avp(:, 4:6) - avpRef(:, 4:6);
    posErr = avp(:, 7:9) - avpRef(:, 7:9);
    posErr(:, 1:2) = posErr(:, 1:2) * glv.Re;
    res.att = sqrt(sum(attErrDeg .^ 2, 2));
    res.vel = sqrt(sum(velErr .^ 2, 2));
    res.pos = sqrt(sum(posErr .^ 2, 2));
    res.t = t;
    detRaw = res.att > detCfg.attDeg | res.vel > detCfg.vel | ...
        res.pos > detCfg.posM;
    detMask = filter(ones(detCfg.minPersist, 1), 1, double(detRaw)) >= ...
        detCfg.minPersist;
end

function metrics = calcFdiMetrics(faultMask, detMask, t)
%CALCFDIMETRICS Calculate FDR, FAR, and detection delay.
%   metrics = CALCFDIMETRICS(faultMask, detMask, t) returns metrics.
    faultMask = logical(faultMask);
    detMask = logical(detMask);
    faultCount = sum(faultMask);
    safeCount = sum(~faultMask);
    metrics.fdr = safeDivide(sum(detMask & faultMask), faultCount);
    metrics.far = safeDivide(sum(detMask & ~faultMask), safeCount);
    metrics.delay = calcDetectionDelay(faultMask, detMask, t);
end

function delay = calcDetectionDelay(faultMask, detMask, t)
%CALCDETECTIONDELAY Mean delay between fault start and detection.
%   delay = CALCDETECTIONDELAY(faultMask, detMask, t) returns mean delay.
    starts = find(diff([0; faultMask]) == 1);
    if isempty(starts)
        delay = 0;
        return;
    end
    delays = nan(numel(starts), 1);
    for k = 1:numel(starts)
        idx = starts(k):numel(faultMask);
        detIdx = find(detMask(idx), 1, 'first');
        if ~isempty(detIdx)
            delays(k) = t(idx(detIdx)) - t(starts(k));
        end
    end
    delay = mean(delays(~isnan(delays)));
end

function out = safeDivide(num, den)
%SAFEDIVIDE Safe divide helper for metric ratios.
%   out = SAFEDIVIDE(num, den) returns num/den or 0 when den is 0.
    if den <= 0
        out = 0;
    else
        out = num / den;
    end
end

function rmse = computeRmse(avp, avpRef)
%COMPUTERMSE Compute RMSE for attitude, velocity, and position.
%   rmse = COMPUTERMSE(avp, avpRef) returns fields att, vel, pos.
    global glv
    attErrDeg = (avp(:, 1:3) - avpRef(:, 1:3)) / glv.deg;
    velErr = avp(:, 4:6) - avpRef(:, 4:6);
    posErr = avp(:, 7:9) - avpRef(:, 7:9);
    posErr(:, 1:2) = posErr(:, 1:2) * glv.Re;
    rmse.att = sqrt(mean(attErrDeg .^ 2, 1));
    rmse.vel = sqrt(mean(velErr .^ 2, 1));
    rmse.pos = sqrt(mean(posErr .^ 2, 1));
end

function plotFdiResults(t, res, detMask, faultMask, detCfg, figTitle)
%PLOTFDIRESULTS Plot residuals and detection indicators.
%   PLOTFDIRESULTS(t, res, detMask, faultMask, detCfg, figTitle) draws curves.
    figure('Name', figTitle);
    subplot(4, 1, 1);
    plot(t, res.att, 'b'); hold on;
    yline(detCfg.attDeg, 'r--');
    plotFaultPatch(t, faultMask);
    ylabel('Att err (deg)');
    title(figTitle);
    subplot(4, 1, 2);
    plot(t, res.vel, 'b'); hold on;
    yline(detCfg.vel, 'r--');
    plotFaultPatch(t, faultMask);
    ylabel('Vel err (m/s)');
    subplot(4, 1, 3);
    plot(t, res.pos, 'b'); hold on;
    yline(detCfg.posM, 'r--');
    plotFaultPatch(t, faultMask);
    ylabel('Pos err (m)');
    subplot(4, 1, 4);
    plot(t, double(detMask), 'k'); hold on;
    plotFaultPatch(t, faultMask);
    xlabel('Time (s)');
    ylabel('Detect');
    ylim([-0.1 1.1]);
end

function plotNavErrors(t, avpNom, avpFault, avpMit, avpRef)
%PLOTNAVERRORS Plot navigation errors for nominal/fault/mitigated cases.
%   PLOTNAVERRORS(t, avpNom, avpFault, avpMit, avpRef) draws curves.
    global glv
    [attNom, velNom, posNom] = navErrors(avpNom, avpRef, glv);
    [attFault, velFault, posFault] = navErrors(avpFault, avpRef, glv);
    [attMit, velMit, posMit] = navErrors(avpMit, avpRef, glv);
    figure('Name', 'Navigation error comparison');
    subplot(3, 1, 1);
    plot(t, attNom(:, 1), 'b', t, attFault(:, 1), 'r', ...
        t, attMit(:, 1), 'g');
    ylabel('Att err (deg)');
    legend('Nominal', 'Faulty', 'Mitigated');
    subplot(3, 1, 2);
    plot(t, velNom(:, 1), 'b', t, velFault(:, 1), 'r', ...
        t, velMit(:, 1), 'g');
    ylabel('Vel err (m/s)');
    subplot(3, 1, 3);
    plot(t, posNom(:, 1), 'b', t, posFault(:, 1), 'r', ...
        t, posMit(:, 1), 'g');
    xlabel('Time (s)');
    ylabel('Pos err (m)');
end

function [attErr, velErr, posErr] = navErrors(avp, avpRef, glv)
%NAVERRORS Compute component-wise navigation errors.
%   [attErr, velErr, posErr] = NAVERRORS(avp, avpRef, glv) returns errors.
    attErr = (avp(:, 1:3) - avpRef(:, 1:3)) / glv.deg;
    velErr = avp(:, 4:6) - avpRef(:, 4:6);
    posErr = avp(:, 7:9) - avpRef(:, 7:9);
    posErr(:, 1:2) = posErr(:, 1:2) * glv.Re;
end

function plotFaultPatch(t, faultMask)
%PLOTFAULTPATCH Add shaded fault window.
%   PLOTFAULTPATCH(t, faultMask) shades fault intervals on plots.
    yl = ylim;
    idx = find(diff([0; faultMask; 0]) ~= 0);
    if numel(idx) < 2
        return;
    end
    for k = 1:2:numel(idx)
        x1 = t(idx(k));
        x2 = t(idx(k + 1) - 1);
        patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end
    uistack(findobj(gca, 'Type', 'line'), 'top');
end
