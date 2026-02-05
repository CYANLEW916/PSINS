% SINS sensor data simulation for IRS1/2 and ISIS based on test_SINS.m.
% Requires 'trj_NLG_approach.mat' from NLG approach trajectory data.
% See also  test_SINS, test_SINS_trj, imuerrset, imuadderr, inspure.

glvs
close all

trjPath = 'D:\All Model\psins251010\data\trj_NLG_approach.mat';
trj = trjfile(trjPath);

%% Initial AVP error setting
initAvpErr = avperrset([0.5; 0.5; 5], 0.1, [10; 10; 10]);
avp00 = avpadderr(trj.avp0, initAvpErr);
trj = bhsimu(trj, 1, 10, 3, trj.ts);

%% Generate "GPS" from reference trajectory (simulating real GPS measurements)
% Downsample reference trajectory to 1Hz GPS rate


gpsRate = 1;  % Hz
gpsInterval = round(1 / (gpsRate * trj.ts));
gpsIdx = 1:gpsInterval:size(trj.avp, 1);

% GPS measurement: [vn(3), pos(3), t] with noise
gpsNoise = struct();
gpsNoise.vel = 0.1;    % m/s (1-sigma)
gpsNoise.pos = 5.0;    % m (1-sigma)

gpsRef = generateGpsMeasurement(trj.avp, gpsIdx, gpsNoise, glv);
fprintf('Generated GPS measurements: %d points at %d Hz\n', ...
    size(gpsRef, 1), gpsRate);

%% Sensor specs
specs = makeSensorSpecs();

%% Simulation with closed-loop correction
sensorData = repmat(struct('name', '', 'imu', [], 'avp', [], 'xkpk', []), ...
    numel(specs), 1);

for k = 1:numel(specs)
    fprintf('\nProcessing %s (closed-loop)...\n', specs(k).name);
    rng(100 + k, 'twister');

    % Generate IMU with errors
    imuerr = imuerrset(specs(k).gyroBias, specs(k).accBias, ...
        specs(k).gyroArw, 0);
    imu = imuadderr(trj.imu, imuerr);
    % Initialize INS
    ins = insinit(avp00, trj.ts);

    % Setup Kalman filter parameters based on sensor grade
    [davp, kfImuerr, rk, fbstr] = setupKfParams(specs(k), glv);

    % Run closed-loop SINS/GPS integration
    [avp, xkpk] = sinsgpsClosedLoop(imu, gpsRef, ins, davp, kfImuerr, ...
        rk, fbstr);

    % Add output noise if specified
    avp = addOutputNoise(avp, specs(k).attNoise, specs(k).posNoise, glv);

    sensorData(k).name = specs(k).name;
    sensorData(k).imu = imu;
    sensorData(k).avp = avp;
    sensorData(k).xkpk = xkpk;  % Store KF state for analysis

    fprintf('%s completed: %d samples\n', specs(k).name, size(avp, 1));
end

%% Also generate open-loop (pure inertial) data for comparison
fprintf('\nGenerating open-loop (pure INS) data for comparison...\n');
sensorDataOpenLoop = repmat(struct('name', '', 'imu', [], 'avp', []), ...
    numel(specs), 1);

for k = 1:numel(specs)
    rng(100 + k, 'twister');
    imuerr = imuerrset(specs(k).gyroBias, specs(k).accBias, ...
        specs(k).gyroArw, 0);
    imu = imuadderr(trj.imu, imuerr);
    avp = inspure(imu, avp00, trj.bh, 0);
    avp = addOutputNoise(avp, specs(k).attNoise, specs(k).posNoise, glv);

    sensorDataOpenLoop(k).name = specs(k).name;
    sensorDataOpenLoop(k).imu = imu;
    sensorDataOpenLoop(k).avp = avp;
end

%% Save results
save('trj10ms_sensor_data.mat', 'sensorData', 'specs', 'trj', 'gpsRef');
save('trj10ms_sensor_data_openloop.mat', ...
    'sensorDataOpenLoop', 'specs', 'trj');
fprintf('\nData saved to trj10ms_sensor_data.mat (closed-loop)\n');
fprintf('Open-loop data saved to trj10ms_sensor_data_openloop.mat\n');

%% Visualization: Compare closed-loop vs open-loop
plotClosedLoopComparison(sensorData, sensorDataOpenLoop, trj, glv);

%% ===================== Helper Functions =====================
function specs = makeSensorSpecs()
%MAKESENSORSPECS Build sensor spec table for IRS1/2 and ISIS.
%   Output:
%     specs - struct array of sensor parameters.
%   IRS1/2: High-precision inertial navigation grade.
%   ISIS:   Lower precision but stable (navigation grade).
    specs = struct('name', {'IRS1', 'IRS2', 'ISIS'}, ...
        'gyroBias', {0.01, 0.01, 0.1}, ...   % deg/h
        'gyroArw', {0.003, 0.003, 0.1}, ...  % deg/sqrt(h)
        'accBias', {100, 100, 5000}, ...     % ug
        'attNoise', {0, 0, 0}, ...           % deg (output noise)
        'posNoise', {0, 0, 0});              % m (output noise)
end

function gps = generateGpsMeasurement(avpRef, idx, noise, glv)
%GENERATEGPSMEASUREMENT Generate GPS measurements from reference trajectory.
%   Inputs:
%     avpRef - reference AVP data.
%     idx    - indices for downsampled GPS epochs.
%     noise  - struct with vel/pos 1-sigma noise settings.
%     glv    - global constants.
%   Output:
%     gps - measurements as [vn(3), pos(3), t].
    nGps = length(idx);
    gps = zeros(nGps, 7);

    for i = 1:nGps
        k = idx(i);
        % Velocity with noise
        vn = avpRef(k, 4:6)' + randn(3, 1) * noise.vel;
        % Position with noise (convert m to rad for lat/lon)
        posNoise = randn(3, 1) * noise.pos;
        posNoise(1:2) = posNoise(1:2) / glv.Re;
        pos = avpRef(k, 7:9)' + posNoise;
        t = avpRef(k, 10);

        gps(i, :) = [vn', pos', t];
    end
end

function [davp, imuerr, rk, fbstr] = setupKfParams(spec, glv)
%SETUPKFPARAMS Setup Kalman filter parameters based on sensor grade.
%   Inputs:
%     spec - sensor specification struct.
%     glv  - global constants.
%   Outputs:
%     davp  - initial AVP error estimates.
%     imuerr - IMU error model.
%     rk    - GPS measurement noise.
%     fbstr - feedback switch string.
    %#ok<NASGU>
    % Initial AVP error estimates (for P0)
    if contains(spec.name, 'ISIS')
        % ISIS: larger initial errors due to lower precision
        davp = avperrset([30; 60], 0.5, [20; 50]);
        imuerr = imuerrset(0.1, 5000, 0.1, 50);
        rk = vperrset(0.2, 10);
        fbstr = 'avped';
    else
        % IRS1/2: smaller initial errors, higher precision
        davp = avperrset([10; 30], 0.3, [15; 30]);
        imuerr = imuerrset(0.01, 100, 0.003, 10);
        rk = vperrset(0.1, 5);
        fbstr = 'avped';
    end
end

function [avp, xkpk] = sinsgpsClosedLoop(imu, gps, ins, davp, imuerr, ...
    rk, fbstr)
%SINSGPSCLOSEDLOOP Simplified closed-loop SINS/GPS integration.
%   Inputs:
%     imu    - IMU data.
%     gps    - GPS measurement data.
%     ins    - initialized INS structure.
%     davp   - initial AVP error estimates.
%     imuerr - IMU error model.
%     rk     - GPS measurement noise.
%     fbstr  - feedback switch string.
%   Outputs:
%     avp  - navigation states over time.
%     xkpk - KF state/history logs.
%
%   This is a simplified version of sinsgps for this demo.
%   Uses 15-state Kalman filter: [phi(3); dvn(3); dpos(3); eb(3); db(3)].

    [nn, ts, nts] = nnts(2, diff(imu(1:2, end)));
    len = size(imu, 1);

    % Initialize Kalman filter (15 states)
    kf = initKf15(davp, imuerr, rk, nts, fbstr);

    % Preallocate output
    nOut = fix(len / nn);
    avp = zeros(nOut, 10);
    xkpk = zeros(nOut, 31);  % [xk(15); diag(Pk)(15); t]

    % GPS index management
    gpsIdx = 1;
    gpsTime = gps(:, 7);

    ki = 1;
    for k = 1:nn:len-nn+1
        k1 = k + nn - 1;
        wvm = imu(k:k1, 1:6);
        t = imu(k1, end);

        % INS mechanization
        ins = insupdate(ins, wvm);

        % KF time update
        kf.Phikk_1 = kffk15(ins, nts);
        kf = kfupdate(kf);

        % Check for GPS measurement
        if gpsIdx <= size(gps, 1) && t >= gpsTime(gpsIdx)
            % GPS measurement available
            gpsVn = gps(gpsIdx, 1:3)';
            gpsPos = gps(gpsIdx, 4:6)';

            % Measurement residual
            zk = [ins.vn - gpsVn; ins.pos - gpsPos];

            % Measurement matrix H
            kf.Hk = [zeros(3), eye(3), zeros(3, 9); ...
                     zeros(3, 6), eye(3), zeros(3, 6)];

            % KF measurement update
            kf = kfupdate(kf, zk, 'M');
            gpsIdx = gpsIdx + 1;
        end

        % Feedback correction
        [kf, ins] = kffeedback15(kf, ins, fbstr, nts);

        % Store results
        avp(ki, :) = [ins.att', ins.vn', ins.pos', t];
        xkpk(ki, :) = [kf.xk', diag(kf.Pxk)', t];
        ki = ki + 1;
    end

    avp(ki:end, :) = [];
    xkpk(ki:end, :) = [];
end

function kf = initKf15(davp, imuerr, rk, nts, fbstr)
%INITKF15 Initialize 15-state Kalman filter.
%   Inputs:
%     davp  - initial AVP error estimates.
%     imuerr - IMU error model.
%     rk    - GPS measurement noise.
%     nts   - INS update interval.
%     fbstr - feedback switch string.
%   Output:
%     kf - initialized Kalman filter structure.
    kf = struct();
    kf.n = 15;
    kf.m = 6;

    % Process noise
    kf.Qt = diag([imuerr.web; imuerr.wdb; zeros(9, 1)]).^2;

    % Measurement noise
    kf.Rk = diag(rk).^2;

    % Initial state covariance
    kf.Pxk = diag([davp; imuerr.eb; imuerr.db]).^2;

    % State and measurement matrices
    kf.xk = zeros(15, 1);
    kf.Hk = zeros(6, 15);

    % Continuous to discrete conversion
    kf.Qk = kf.Qt * nts;
    kf.Gammak = eye(15);
    kf.fbstr = fbstr;
    kf = kfinit0(kf, nts);
end

function Phi = kffk15(ins, nts)
%KFFK15 Compute state transition matrix for 15-state filter.
%   Inputs:
%     ins - INS structure.
%     nts - time step.
%   Output:
%     Phi - discrete-time transition matrix.
    eth = ins.eth;
    Cnb = ins.Cnb;
    fn = ins.fn;

    % Simplified state transition (first-order approximation)
    Maa = askew(-eth.wnin);
    Mav = [0, -1 / eth.RMh, 0; 1 / eth.clRNh, 0, 0; ...
        tan(ins.pos(1)) / eth.clRNh, 0, 0];
    Map = zeros(3);
    Mva = askew(fn);
    Mvv = askew(-2 * eth.wnie - eth.wnen);
    Mvp = zeros(3);
    Mpv = [0, 1 / eth.RMh, 0; 1 / eth.clRNh, 0, 0; 0, 0, 1];
    Mpp = zeros(3);

    % State transition matrix
    F = [Maa, Mav, Map, -Cnb, zeros(3); ...
         Mva, Mvv, Mvp, zeros(3), Cnb; ...
         zeros(3), Mpv, Mpp, zeros(3, 6); ...
         zeros(6, 15)];

    Phi = eye(15) + F * nts;
end

function [kf, ins] = kffeedback15(kf, ins, fbstr, nts)
%KFFEEDBACK15 Feedback correction for 15-state filter.
%   Inputs:
%     kf    - Kalman filter structure.
%     ins   - INS structure.
%     fbstr - feedback switch string.
%     nts   - time step.
%   Outputs:
%     kf  - updated Kalman filter structure.
%     ins - corrected INS structure.
    %#ok<INUSD>
    % Extract state corrections
    phi = kf.xk(1:3);
    dvn = kf.xk(4:6);
    dpos = kf.xk(7:9);
    deb = kf.xk(10:12);
    ddb = kf.xk(13:15);

    % Feedback gain (exponential decay)
    beta = 1.0;  % Full feedback

    % Apply corrections
    if contains(fbstr, 'a')
        ins.qnb = qdelphi(ins.qnb, beta * phi);
        [ins.qnb, ins.att, ins.Cnb] = attsyn(ins.qnb);
        kf.xk(1:3) = (1 - beta) * phi;
    end

    if contains(fbstr, 'v')
        ins.vn = ins.vn - beta * dvn;
        kf.xk(4:6) = (1 - beta) * dvn;
    end

    if contains(fbstr, 'p')
        ins.pos = ins.pos - beta * dpos;
        kf.xk(7:9) = (1 - beta) * dpos;
    end

    if contains(fbstr, 'e')
        ins.eb = ins.eb + beta * deb;
        kf.xk(10:12) = (1 - beta) * deb;
    end

    if contains(fbstr, 'd')
        ins.db = ins.db + beta * ddb;
        kf.xk(13:15) = (1 - beta) * ddb;
    end
end

function avp = addOutputNoise(avp, attNoiseDeg, posNoiseM, glv)
%ADDOUTPUTNOISE Add attitude/position output noise to AVP.
%   Inputs:
%     avp         - AVP output array.
%     attNoiseDeg - attitude noise in degrees.
%     posNoiseM   - position noise in meters.
%     glv         - global constants.
%   Output:
%     avp - noisy AVP output.
    if attNoiseDeg > 0
        attNoiseRad = attNoiseDeg * glv.deg;
        attPert = randn(size(avp, 1), 3) * attNoiseRad;
        avp(:, 1:3) = avp(:, 1:3) + attPert;
    end
    if posNoiseM > 0
        posPert = randn(size(avp, 1), 3) * posNoiseM;
        posPert(:, 1:2) = posPert(:, 1:2) / glv.Re;
        avp(:, 7:9) = avp(:, 7:9) + posPert;
    end
end

function plotClosedLoopComparison(sensorCL, sensorOL, trj, glv)
%PLOTCLOSEDLOOPCOMPARISON Compare closed-loop vs open-loop navigation.
%   Inputs:
%     sensorCL - closed-loop sensor data.
%     sensorOL - open-loop sensor data.
%     trj      - reference trajectory.
%     glv      - global constants.
    figure('Name', 'Closed-Loop vs Open-Loop Comparison', ...
        'Position', [50, 50, 1400, 900]);

    avpRef = trj.avp;
    tRef = avpRef(:, 10) - avpRef(1, 10);

    colors = {'b', 'r', 'm'};
    sensorNames = {'IRS1', 'IRS2', 'ISIS'};

    % Position error comparison
    for i = 1:3
        subplot(3, 3, i);
        hold on;

        % Closed-loop error
        avpCL = sensorCL(i).avp;
        tCL = avpCL(:, 10) - avpCL(1, 10);
        errCL = computePosError(avpCL, avpRef, glv);

        % Open-loop error
        avpOL = sensorOL(i).avp;
        tOL = avpOL(:, 10) - avpOL(1, 10);
        errOL = computePosError(avpOL, avpRef, glv);

        plot(tOL, errOL, [colors{i}, '--'], 'LineWidth', 1);
        plot(tCL, errCL, [colors{i}, '-'], 'LineWidth', 1.5);

        xlabel('Time (s)');
        ylabel('Position Error (m)');
        title(sprintf('%s Position Error', sensorNames{i}));
        legend({'Open-Loop', 'Closed-Loop'}, 'Location', 'best');
        grid on;
    end

    % Attitude error comparison
    for i = 1:3
        subplot(3, 3, 3 + i);
        hold on;

        avpCL = sensorCL(i).avp;
        tCL = avpCL(:, 10) - avpCL(1, 10);
        errCL = computeAttError(avpCL, avpRef, glv);

        avpOL = sensorOL(i).avp;
        tOL = avpOL(:, 10) - avpOL(1, 10);
        errOL = computeAttError(avpOL, avpRef, glv);

        plot(tOL, errOL, [colors{i}, '--'], 'LineWidth', 1);
        plot(tCL, errCL, [colors{i}, '-'], 'LineWidth', 1.5);

        xlabel('Time (s)');
        ylabel('Attitude Error (deg)');
        title(sprintf('%s Attitude Error', sensorNames{i}));
        legend({'Open-Loop', 'Closed-Loop'}, 'Location', 'best');
        grid on;
    end

    % Trajectory comparison
    subplot(3, 3, [7, 8, 9]);
    hold on;

    lat0 = avpRef(1, 7);
    lon0 = avpRef(1, 8);

    % Reference
    xRef = (avpRef(:, 8) - lon0) * glv.Re * cos(lat0);
    yRef = (avpRef(:, 7) - lat0) * glv.Re;
    plot(xRef, yRef, 'k-', 'LineWidth', 2.5);

    % Closed-loop sensors
    for i = 1:3
        avpCL = sensorCL(i).avp;
        xCL = (avpCL(:, 8) - lon0) * glv.Re * cos(lat0);
        yCL = (avpCL(:, 7) - lat0) * glv.Re;
        plot(xCL, yCL, [colors{i}, '-'], 'LineWidth', 1.5);
    end

    % Open-loop ISIS (to show drift)
    avpOL = sensorOL(3).avp;
    xOL = (avpOL(:, 8) - lon0) * glv.Re * cos(lat0);
    yOL = (avpOL(:, 7) - lat0) * glv.Re;
    plot(xOL, yOL, 'm--', 'LineWidth', 1);

    xlabel('East (m)');
    ylabel('North (m)');
    title('Trajectory Comparison');
    legend({'Reference', 'IRS1 (CL)', 'IRS2 (CL)', 'ISIS (CL)', ...
        'ISIS (OL)'}, 'Location', 'best');
    axis equal;
    grid on;

    % Print summary statistics
    fprintf('\n========== Position RMSE Summary ==========');
    fprintf('\nSensor      Open-Loop(m)   Closed-Loop(m)   Improvement');
    for i = 1:3
        errOL = computePosError(sensorOL(i).avp, avpRef, glv);
        errCL = computePosError(sensorCL(i).avp, avpRef, glv);
        rmseOL = sqrt(mean(errOL .^ 2));
        rmseCL = sqrt(mean(errCL .^ 2));
        improvement = (1 - rmseCL / rmseOL) * 100;
        fprintf('\n%-8s    %10.1f     %10.1f       %5.1f%%', ...
            sensorNames{i}, rmseOL, rmseCL, improvement);
    end
    fprintf('\n');
end

function err = computePosError(avp, avpRef, glv)
%COMPUTEPOSERROR Compute position error magnitude in meters.
%   Inputs:
%     avp    - AVP data to evaluate.
%     avpRef - reference AVP data.
%     glv    - global constants.
%   Output:
%     err - position error magnitude in meters.
    tAvp = avp(:, 10);
    tRef = avpRef(:, 10);

    % Interpolate reference to match AVP times
    posRef = interp1(tRef, avpRef(:, 7:9), tAvp, 'linear', 'extrap');
    posAvp = avp(:, 7:9);

    % Convert to meters
    lat0 = mean(posRef(:, 1));
    dpos = posAvp - posRef;
    dposM = [dpos(:, 1) * glv.Re, dpos(:, 2) * glv.Re * cos(lat0), ...
        dpos(:, 3)];

    err = sqrt(sum(dposM .^ 2, 2));
end

function err = computeAttError(avp, avpRef, glv)
%COMPUTEATTERROR Compute attitude error magnitude in degrees.
%   Inputs:
%     avp    - AVP data to evaluate.
%     avpRef - reference AVP data.
%     glv    - global constants.
%   Output:
%     err - attitude error magnitude in degrees.
    tAvp = avp(:, 10);
    tRef = avpRef(:, 10);

    attRef = interp1(tRef, avpRef(:, 1:3), tAvp, 'linear', 'extrap');
    attAvp = avp(:, 1:3);

    datt = (attAvp - attRef) / glv.deg;
    datt(:, 3) = mod(datt(:, 3) + 180, 360) - 180;  % Wrap yaw

    err = sqrt(sum(datt .^ 2, 2));
end
