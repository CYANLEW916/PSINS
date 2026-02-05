% Fault injection and detection evaluation for IRS1/2 sensors.
% Injects detectable step faults and compares FDR/FAR/delay and RMSE results
% with/without detection mitigation.
% Requires 'trj10ms_sensor_data.mat' from test_SINS_IRS_ISIS.m before running.
% See also  test_SINS_IRS_ISIS, test_SINS_trj.

glvs

if exist('trj10ms_sensor_data.mat', 'file') ~= 2
    test_SINS_IRS_ISIS;
end

load('trj10ms_sensor_data.mat', 'sensorData', 'specs', 'trj');

irsNames = {'IRS1', 'IRS2'};
irsIdx = zeros(numel(irsNames), 1);
for k = 1:numel(irsNames)
    irsIdx(k) = find(strcmp({sensorData.name}, irsNames{k}), 1, 'first');
    if isempty(irsIdx(k))
        error('IRS sensor data not found: %s.', irsNames{k});
    end
end

imuNom = sensorData(irsIdx(1)).imu;
t = imuNom(:, end);
nSample = size(imuNom, 1);

faultCfg = buildFaultConfig(t);
detCfg = makeDetectionConfig();

[refImuAll, resNomAll] = computeAllImuResiduals(sensorData, nSample, t);
detCfg = tuneDetectionThresholds(resNomAll, detCfg);
results = struct([]);
for k = 1:numel(irsIdx)
    idx = irsIdx(k);
    imuNom = sensorData(idx).imu;
    [imuFault, faultMask, statusMask] = injectImuFaults(imuNom, t, faultCfg);
    sensorDataFault = sensorData;
    sensorDataFault(idx).imu = imuFault;

    [resFault, detMask] = detectImuFaults(imuFault, refImuAll{idx}, t, detCfg);
    [resNom, detMaskNom] = detectImuFaults(imuNom, refImuAll{idx}, t, detCfg);
    metricsFault = calcFdiMetrics(faultMask, detMask, t);
    metricsNom = calcFdiMetrics(false(size(t)), detMaskNom, t);

    avpNom = sensorData(idx).avp;
    avpCount = size(avpNom, 1);
    tAvp = t(1:avpCount);
    avpFault = inspure(imuFault, avpNom(1, 1:9)', trj.bh, 1);
    avpFault = avpFault(1:avpCount, :);
    sensorDataFault(idx).avp = avpFault;

    imuMit = imuFault;
    imuMit(detMask, :) = refImuAll{idx}(detMask, :);
    avpMit = inspure(imuMit, avpNom(1, 1:9)', trj.bh, 1);
    avpMit = avpMit(1:avpCount, :);

    trjAvp = trj.avp(1:avpCount, :);
    rmseNom = computeRmse(avpNom, trjAvp);
    rmseFault = computeRmse(avpFault, trjAvp);
    rmseMit = computeRmse(avpMit, trjAvp);

    fprintf(['\n[%s] Fault detection metrics (faulty): FDR=%.3f, ', ...
        'FAR=%.3f, delay=%.2fs\n'], sensorData(idx).name, metricsFault.fdr, ...
        metricsFault.far, metricsFault.delay);
    fprintf(['[%s] Fault detection metrics (nominal): FDR=%.3f, ', ...
        'FAR=%.3f, delay=%.2fs\n'], sensorData(idx).name, metricsNom.fdr, ...
        metricsNom.far, metricsNom.delay);
    fprintf(['[%s] RMSE att/vel/pos (nominal): [%.3f %.3f %.3f] ', ...
        '[%.3f %.3f %.3f] [%.3f %.3f %.3f]\n'], sensorData(idx).name, ...
        rmseNom.att, rmseNom.vel, rmseNom.pos);
    fprintf(['[%s] RMSE att/vel/pos (faulty): [%.3f %.3f %.3f] ', ...
        '[%.3f %.3f %.3f] [%.3f %.3f %.3f]\n'], sensorData(idx).name, ...
        rmseFault.att, rmseFault.vel, rmseFault.pos);
    fprintf(['[%s] RMSE att/vel/pos (mitigated): [%.3f %.3f %.3f] ', ...
        '[%.3f %.3f %.3f] [%.3f %.3f %.3f]\n'], sensorData(idx).name, ...
        rmseMit.att, rmseMit.vel, rmseMit.pos);

    results(k).name = sensorData(idx).name;
    results(k).faultMask = faultMask;
    results(k).statusMask = statusMask;
    results(k).metricsFault = metricsFault;
    results(k).metricsNom = metricsNom;
    results(k).rmseNom = rmseNom;
    results(k).rmseFault = rmseFault;
    results(k).rmseMit = rmseMit;
    results(k).avpNom = avpNom;
    results(k).avpFault = avpFault;
    results(k).avpMit = avpMit;
    results(k).tAvp = tAvp;
    results(k).resFault = resFault;
    results(k).resNom = resNom;
    results(k).detMask = detMask;
    results(k).detMaskNom = detMaskNom;

    plotFdiResults(t, resFault, detMask, faultMask, detCfg, ...
        sprintf('%s faulty case', sensorData(idx).name));
    plotFdiResults(t, resNom, detMaskNom, false(size(t)), detCfg, ...
        sprintf('%s nominal case', sensorData(idx).name));
    plotNavErrors(tAvp, avpNom, avpFault, avpMit, trjAvp, sensorData(idx).name);
end

save('trj10ms_sensor_data_faults.mat', 'specs', 'faultCfg', 'detCfg', ...
    'results', 't', 'refImuAll');

%% helper functions
function faultCfg = buildFaultConfig(t)
%BUILDFAULTCONFIG Configure detectable faults for IRS sensors.
%   faultCfg = BUILDFAULTCONFIG(t) returns a struct array with fault windows
%   and magnitudes based on the time vector t.
    tEnd = t(end);
    faultCfg = struct( ...
        'name', {'gyro_step', 'acc_step'}, ...
        'tStart', {0.2 * tEnd, 0.6 * tEnd}, ...
        'tEnd', {0.35 * tEnd, 0.75 * tEnd}, ...
        'gyroDeg', {[0.6 -0.5 0.4], []}, ...
        'acc', {[], [0.03 -0.02 0.015]});
end

function refImu = buildReferenceImu(sensorData, faultIdx, nSample)
%BUILDREFERENCEIMU Build reference IMU from non-faulty sensors.
%   refImu = BUILDREFERENCEIMU(sensorData, faultIdx, nSample) returns a
%   reference IMU using the mean of the other sensors.
    otherIdx = setdiff(1:numel(sensorData), faultIdx);
    if isempty(otherIdx)
        error('No reference sensor data available.');
    end
    imuStack = [];
    for k = otherIdx
        imu = sensorData(k).imu(1:nSample, :);
        imuStack = cat(3, imuStack, imu);
    end
    refImu = mean(imuStack, 3);
end

function [refImuAll, resNomAll] = computeAllImuResiduals(sensorData, nSample, t)
%COMPUTEALLIMURESIDUALS Build references and residuals for all sensors.
%   [refImuAll, resNomAll] = COMPUTEALLIMURESIDUALS(sensorData, nSample, t)
%   returns cell arrays for each sensor.
    nSensor = numel(sensorData);
    refImuAll = cell(nSensor, 1);
    resNomAll = cell(nSensor, 1);
    for k = 1:nSensor
        refImuAll{k} = buildReferenceImu(sensorData, k, nSample);
        imu = sensorData(k).imu(1:nSample, :);
        resNomAll{k} = computeImuResiduals(imu, refImuAll{k}, t);
    end
end

function detCfg = makeDetectionConfig()
%MAKEDETECTIONCONFIG Thresholds for fault detection.
%   detCfg = MAKEDETECTIONCONFIG() returns thresholds and persistence.
    detCfg = struct('sigmaFactor', 3, 'baseFrac', 0.15, ...
        'minPersist', 5, 'gyroDeg', [], 'acc', []);
end

function [imuFault, faultMask, statusMask] = injectImuFaults(imu, t, faultCfg)
%INJECTIMUFAULTS Inject additive faults into IMU with recovery.
%   [imuFault, faultMask, statusMask] = INJECTIMUFAULTS(imu, t, faultCfg)
%   injects faults defined in faultCfg into IMU, resets to nominal after
%   each fault window, and returns the fault/status masks.
    global glv
    imuNom = imu;
    imuFault = imuNom;
    faultMask = false(size(t));
    statusMask = false(size(t));
    for k = 1:numel(faultCfg)
        idx = t >= faultCfg(k).tStart & t <= faultCfg(k).tEnd;
        faultMask = faultMask | idx;
        statusMask(idx) = true;
        if ~isempty(faultCfg(k).gyroDeg)
            gyroAdd = repmat(faultCfg(k).gyroDeg, sum(idx), 1) * glv.deg;
            imuFault(idx, 1:3) = imuNom(idx, 1:3) + gyroAdd;
        end
        if ~isempty(faultCfg(k).acc)
            accAdd = repmat(faultCfg(k).acc, sum(idx), 1);
            imuFault(idx, 4:6) = imuNom(idx, 4:6) + accAdd;
        end
    end
end

function res = computeImuResiduals(imu, imuRef, t)
%COMPUTEIMURESIDUALS Compute residual magnitudes between IMU and reference.
%   res = COMPUTEIMURESIDUALS(imu, imuRef, t) returns residuals.
    global glv
    gyroErr = (imu(:, 1:3) - imuRef(:, 1:3)) / glv.deg;
    accErr = imu(:, 4:6) - imuRef(:, 4:6);
    res.gyro = sqrt(sum(gyroErr .^ 2, 2));
    res.acc = sqrt(sum(accErr .^ 2, 2));
    res.t = t;
end

function detCfg = tuneDetectionThresholds(resNom, detCfg)
%TUNEDETECTIONTHRESHOLDS Estimate thresholds from nominal residuals.
%   detCfg = TUNEDETECTIONTHRESHOLDS(resNom, detCfg) updates thresholds.
    t = resNom{1}.t;
    baseIdx = t <= detCfg.baseFrac * t(end);
    if ~any(baseIdx)
        baseIdx = true(size(t));
    end
    gyroStack = [];
    accStack = [];
    for k = 1:numel(resNom)
        gyroStack = [gyroStack; resNom{k}.gyro(baseIdx)];
        accStack = [accStack; resNom{k}.acc(baseIdx)];
    end
    detCfg.gyroDeg = mean(gyroStack) + detCfg.sigmaFactor * std(gyroStack);
    detCfg.acc = mean(accStack) + detCfg.sigmaFactor * std(accStack);
end

function [res, detMask] = detectImuFaults(imu, imuRef, t, detCfg)
%DETECTIMUFAULTS Detect faults based on IMU residual thresholds.
%   [res, detMask] = DETECTIMUFAULTS(imu, imuRef, t, detCfg) returns residuals
%   and the detection mask.
    res = computeImuResiduals(imu, imuRef, t);
    detRaw = res.gyro > detCfg.gyroDeg | res.acc > detCfg.acc;
    detMask = filter(ones(detCfg.minPersist, 1), 1, double(detRaw)) >= ...
        detCfg.minPersist;
end

function [voteMask, voteDetail] = voteFaults(resNomAll, detCfg)
%VOTEFAULTS Vote-based fault flags across sensors.
%   [voteMask, voteDetail] = VOTEFAULTS(resNomAll, detCfg) flags a sensor
%   as faulty when it disagrees with its reference while others agree.
    nSensor = numel(resNomAll);
    voteMask = false(numel(resNomAll{1}.t), nSensor);
    voteDetail = struct('gyro', [], 'acc', []);
    voteDetail.gyro = false(size(voteMask));
    voteDetail.acc = false(size(voteMask));
    for k = 1:nSensor
        res = resNomAll{k};
        voteDetail.gyro(:, k) = res.gyro > detCfg.gyroDeg;
        voteDetail.acc(:, k) = res.acc > detCfg.acc;
        voteMask(:, k) = voteDetail.gyro(:, k) | voteDetail.acc(:, k);
    end
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
    subplot(3, 1, 1);
    plot(t, res.gyro, 'b'); hold on;
    yline(detCfg.gyroDeg, 'r--');
    plotFaultPatch(t, faultMask);
    ylabel('Gyro err (deg)');
    title(figTitle);
    subplot(3, 1, 2);
    plot(t, res.acc, 'b'); hold on;
    yline(detCfg.acc, 'r--');
    plotFaultPatch(t, faultMask);
    ylabel('Acc err (m/s)');
    subplot(3, 1, 3);
    plot(t, double(detMask), 'k'); hold on;
    plotFaultPatch(t, faultMask);
    xlabel('Time (s)');
    ylabel('Detect');
    ylim([-0.1 1.1]);
end

function plotNavErrors(t, avpNom, avpFault, avpMit, avpRef, sensorName)
%PLOTNAVERRORS Plot navigation errors for nominal/fault/mitigated cases.
%   PLOTNAVERRORS(t, avpNom, avpFault, avpMit, avpRef, sensorName) draws
%   curves for the specified sensor.
    global glv
    [attNom, velNom, posNom] = navErrors(avpNom, avpRef, glv);
    [attFault, velFault, posFault] = navErrors(avpFault, avpRef, glv);
    [attMit, velMit, posMit] = navErrors(avpMit, avpRef, glv);
    attNomNorm = vecnorm(attNom, 2, 2);
    attFaultNorm = vecnorm(attFault, 2, 2);
    attMitNorm = vecnorm(attMit, 2, 2);
    velNomNorm = vecnorm(velNom, 2, 2);
    velFaultNorm = vecnorm(velFault, 2, 2);
    velMitNorm = vecnorm(velMit, 2, 2);
    posNomNorm = vecnorm(posNom, 2, 2);
    posFaultNorm = vecnorm(posFault, 2, 2);
    posMitNorm = vecnorm(posMit, 2, 2);
    figure('Name', sprintf('Navigation error comparison (%s)', sensorName));
    subplot(3, 1, 1);
    plot(t, attNomNorm, 'b', t, attFaultNorm, 'r', t, attMitNorm, 'g');
    ylabel('Att err (deg)');
    legend('Nominal', 'Faulty', 'Mitigated');
    subplot(3, 1, 2);
    plot(t, velNomNorm, 'b', t, velFaultNorm, 'r', t, velMitNorm, 'g');
    ylabel('Vel err (m/s)');
    subplot(3, 1, 3);
    plot(t, posNomNorm, 'b', t, posFaultNorm, 'r', t, posMitNorm, 'g');
    xlabel('Time (s)');
    ylabel('Pos err (m)');
end

function plotSensorOutputs(t, sensorDataFault, nSample)
%PLOTSENSOROUTPUTS Plot attitude/velocity/position from all sensors.
%   PLOTSENSOROUTPUTS(t, sensorDataFault, nSample) draws three figures with
%   attitude, velocity, and position for the fault-injected sensors.
    global glv
    names = {sensorDataFault.name};
    colors = lines(numel(sensorDataFault));
    figure('Name', 'Attitude outputs (fault injected)');
    subplot(3, 1, 1);
    hold on;
    for k = 1:numel(sensorDataFault)
        avp = sensorDataFault(k).avp(1:nSample, :);
        plot(t, avp(:, 1) / glv.deg, 'Color', colors(k, :));
    end
    ylabel('Pitch (deg)');
    legend(names, 'Location', 'best');
    subplot(3, 1, 2);
    hold on;
    for k = 1:numel(sensorDataFault)
        avp = sensorDataFault(k).avp(1:nSample, :);
        plot(t, avp(:, 2) / glv.deg, 'Color', colors(k, :));
    end
    ylabel('Roll (deg)');
    subplot(3, 1, 3);
    hold on;
    for k = 1:numel(sensorDataFault)
        avp = sensorDataFault(k).avp(1:nSample, :);
        plot(t, avp(:, 3) / glv.deg, 'Color', colors(k, :));
    end
    ylabel('Yaw (deg)');
    xlabel('Time (s)');

    figure('Name', 'Velocity outputs (fault injected)');
    subplot(3, 1, 1);
    hold on;
    for k = 1:numel(sensorDataFault)
        avp = sensorDataFault(k).avp(1:nSample, :);
        plot(t, avp(:, 4), 'Color', colors(k, :));
    end
    ylabel('V_E (m/s)');
    legend(names, 'Location', 'best');
    subplot(3, 1, 2);
    hold on;
    for k = 1:numel(sensorDataFault)
        avp = sensorDataFault(k).avp(1:nSample, :);
        plot(t, avp(:, 5), 'Color', colors(k, :));
    end
    ylabel('V_N (m/s)');
    subplot(3, 1, 3);
    hold on;
    for k = 1:numel(sensorDataFault)
        avp = sensorDataFault(k).avp(1:nSample, :);
        plot(t, avp(:, 6), 'Color', colors(k, :));
    end
    ylabel('V_U (m/s)');
    xlabel('Time (s)');

    figure('Name', 'Position outputs (fault injected)');
    subplot(3, 1, 1);
    hold on;
    for k = 1:numel(sensorDataFault)
        avp = sensorDataFault(k).avp(1:nSample, :);
        plot(t, avp(:, 7) / glv.deg, 'Color', colors(k, :));
    end
    ylabel('Lat (deg)');
    legend(names, 'Location', 'best');
    subplot(3, 1, 2);
    hold on;
    for k = 1:numel(sensorDataFault)
        avp = sensorDataFault(k).avp(1:nSample, :);
        plot(t, avp(:, 8) / glv.deg, 'Color', colors(k, :));
    end
    ylabel('Lon (deg)');
    subplot(3, 1, 3);
    hold on;
    for k = 1:numel(sensorDataFault)
        avp = sensorDataFault(k).avp(1:nSample, :);
        plot(t, avp(:, 9), 'Color', colors(k, :));
    end
    ylabel('Hgt (m)');
    xlabel('Time (s)');
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
