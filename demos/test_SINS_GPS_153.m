% SINS/GPS integrated navigation simulation with residual-based fault detection.
% Please run 'test_SINS_trj.m' to generate 'trj10ms.mat' beforehand!!!
% See also  test_SINS_trj, test_SINS, test_SINS_GPS_186, test_SINS_GPS_193.
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 17/06/2011, updated 15/11/2024
glvs
psinstypedef(153);
trj = trjfile('trj_TL_approach.mat');

% initial settings
[nn, ts, nts] = nnts(2, trj.ts);
imuerr = imuerrset(0.03, 100, 0.001, 5);
davp0 = avperrset([0.5;-0.5;20], 0.1, [1;1;3]);

% fault schedule: [start[s], stop[s], dNorth[m], dEast[m], dUp[m]]
faultWindows = [
    120, 150,  80,  -60,   0;    % sustained horizontal bias
    220, 240,   0,    0, 120;    % vertical bias
];
% impulsive spikes: [time[s], dNorth[m], dEast[m], dUp[m]]
faultSpikes = [
    310, 250, -200,  40;
    420,-180,  220, -35;
];

% Monte Carlo configuration
nMonteCarlo = 30;                      % number of Monte Carlo experiments
rng('shuffle');
simResults = cell(nMonteCarlo,1);

for expIdx = 1:nMonteCarlo
    simResults{expIdx} = runSinsGpsExperiment(trj, imuerr, davp0, faultWindows, faultSpikes, nn, ts, expIdx, nMonteCarlo);
end

firstRun = simResults{1};

% show detailed results for the first Monte Carlo experiment
insplot(firstRun.avpFD);  % best-performing filter
avpcmpplot(trj.avp, firstRun.avpAll, firstRun.avpFD);
avperrFD = avpcmp(firstRun.avpFD, trj.avp);
kfplot(firstRun.xkpkFD, avperrFD, imuerr);

% Position error analysis for the first experiment
fprintf('\nPosition RMSE (All measurements) [North East Up] / m: %6.2f %6.2f %6.2f\n', firstRun.rmseAll);
fprintf('Position RMSE (FD-enabled)      [North East Up] / m: %6.2f %6.2f %6.2f\n', firstRun.rmseFD);
fprintf('3D position RMSE - All: %6.2f m, FD-enabled: %6.2f m\n\n', firstRun.rmseAll3D, firstRun.rmseFD3D);

myfigure('Position Error Comparison');
labels = {'North', 'East', 'Up'};
for idx = 1:3
    subplot(3,1,idx);
    plot(firstRun.tAll, firstRun.posErrAll(:,idx), 'b', ...
         firstRun.tFD, firstRun.posErrFD(:,idx), 'r--', 'LineWidth', 1.2);
    grid on;  ylabel([labels{idx}, ' / m']);
    if idx==1
        title('Position errors with and without residual-based fault rejection (run 1)');
    end
    if idx==3
        xlabel('Time / s');
    end
    legend('All measurements','FD-enabled','Location','best');
end

% residual-based fault detection metrics
nisThreshold = firstRun.nisThreshold;
fdActive = firstRun.detectLog | firstRun.dwellLog;
myfigure('Residual-based Fault Detection');
subplot(2,1,1);
plot(firstRun.tMeas, firstRun.nisLog, 'b', 'LineWidth', 1.2); hold on;
yline(nisThreshold, 'r--', 'LineWidth', 1.2);
grid on;  xlabel('Time / s');  ylabel('NIS');
legend('NIS','Threshold','Location','best');
title('Normalized innovation squared (NIS) test - run 1');
subplot(2,1,2);
stairs(firstRun.tMeas, firstRun.faultFlagLog, 'k', 'LineWidth', 1.2); hold on;
stairs(firstRun.tMeas, fdActive, 'r--', 'LineWidth', 1.2);
stairs(firstRun.tMeas, firstRun.dwellLog, 'Color',[0.3 0.6 0.9], 'LineWidth', 1.1);
grid on;  xlabel('Time / s');
ylabel('State');
ylim([-0.1 1.1]);
legend('Injected fault','FD decision','Dwell active','Location','best');

% Actual Navigation Performance (ANP) estimation versus realised horizontal error (run 1)
myfigure('ANP_vs_HorizontalError');
plot(firstRun.tFD, firstRun.anp, 'b', 'LineWidth', 1.5); hold on;
plot(firstRun.tFD, firstRun.horErrFD, '--', 'Color',[0.1 0.6 0.1], 'LineWidth', 1.5);
grid on;  xygo('Time / s', 'Horizontal performance / m');
legend('ANP (95%)','Horizontal error','Location','best');
title('ANP estimate versus realised horizontal error (run 1)');

% Monte Carlo analysis of ANP confidence
totalSamples = 0;
totalOutOfLimit = 0;
myfigure('ANP_MonteCarlo_Comparison');
hold on;
grid on;
hBgAnp = []; hBgErr = [];
for expIdx = 1:nMonteCarlo
    sim = simResults{expIdx};
    totalSamples = totalSamples + numel(sim.anp);
    totalOutOfLimit = totalOutOfLimit + sum(sim.horErrFD > sim.anp);
    if expIdx == 1
        continue;
    elseif expIdx == 2
        hBgAnp = plot(sim.tFD, sim.anp, 'Color', [0.7 0.82 0.94], 'LineWidth', 0.9);
        hBgErr = plot(sim.tFD, sim.horErrFD, '--', 'Color', [0.76 0.88 0.76], 'LineWidth', 0.9);
    else
        plot(sim.tFD, sim.anp, 'Color', [0.7 0.82 0.94], 'LineWidth', 0.9);
        plot(sim.tFD, sim.horErrFD, '--', 'Color', [0.76 0.88 0.76], 'LineWidth', 0.9);
    end
end
hFirstAnp = plot(firstRun.tFD, firstRun.anp, 'b', 'LineWidth', 1.6);
hFirstErr = plot(firstRun.tFD, firstRun.horErrFD, '--', 'Color',[0.1 0.6 0.1], 'LineWidth', 1.6);
xygo('Time / s', 'Horizontal performance / m');
title('Monte Carlo comparison of ANP against horizontal error');
if nMonteCarlo > 1 && ~isempty(hBgAnp) && ~isempty(hBgErr)
    if nMonteCarlo == 2
        bgLabel = 'run 2';
    else
        bgLabel = sprintf('runs 2-%d', nMonteCarlo);
    end
    legend([hBgAnp, hBgErr, hFirstAnp, hFirstErr], ...
        {['ANP (', bgLabel, ')'], ['Horizontal error (', bgLabel, ')'], ...
         'ANP (run 1)','Horizontal error (run 1)'}, ...
        'Location','best');
else
    legend([hFirstAnp, hFirstErr], {'ANP (run 1)','Horizontal error (run 1)'}, 'Location','best');
end

confidenceLevel = 1 - totalOutOfLimit / max(totalSamples, 1);
fprintf('\nANP Algorithm Monte Carlo Simulation Confidence Level Statistical Table\n');
fprintf('%-24s %-20s %-24s %-12s\n', 'Number of experiments', 'Total number of data', ...
    'Number of out-of-limit cases', 'Confidence');
fprintf('%-24d %-20d %-24d %10.2f%%\n\n', nMonteCarlo, totalSamples, totalOutOfLimit, 100*confidenceLevel);


function sim = runSinsGpsExperiment(trj, imuerr, davp0, faultWindows, faultSpikes, nn, ts, runIdx, nMonteCarlo)
% Execute one Monte Carlo run of the SINS/GPS simulation and return diagnostic data.
    imu = imuadderr(trj.imu, imuerr);
    insAll = insinit(avpadderr(trj.avp0,davp0), ts);
    insFD  = insinit(avpadderr(trj.avp0,davp0), ts);

    % Kalman filters
    rk = poserrset([1;1;3]);
    kfAll = kfinit(insAll, davp0, imuerr, rk);
    kfFD  = kfinit(insFD,  davp0, imuerr, rk);
    commonPmin = [avperrset(0.01,1e-4,0.1); gabias(1e-3, [1,10])].^2;
    kfAll.Pmin = commonPmin;  kfAll.pconstrain = 1;
    kfFD.Pmin  = commonPmin;  kfFD.pconstrain  = 1;

    % residual-based fault detection configuration (NIS gate with 99% chi-square threshold)
    kfFD.fd.enable = 1;
    kfFD.fd.chi2Threshold = 11.345;  % ~chi2inv(0.99, 3)
    kfFD.fd.holdTime = 5;            % dwell time (s) before re-enabling rejected GPS measurements

    len = length(imu);
    rows = fix(len/nn);
    [avpAll, avpFD, xkpkAll, xkpkFD, nisLog, faultFlagLog, detectLog, dwellLog] = ...
        prealloc(rows, 10, 10, 2*kfAll.n+1, 2*kfFD.n+1, 1, 1, 1, 1);
    timebar(nn, len, sprintf('15-state SINS/GPS Simulation (%d/%d).', runIdx, nMonteCarlo));
    ki = 1;

    for k = 1:nn:len-nn+1
        k1 = k+nn-1;
        wvm = imu(k:k1,1:6);  t = imu(k1,end);
        insAll = insupdate(insAll, wvm);
        insFD  = insupdate(insFD,  wvm);
        kfAll.Phikk_1 = kffk(insAll);
        kfFD.Phikk_1  = kffk(insFD);
        kfAll = kfupdate(kfAll);
        kfFD  = kfupdate(kfFD);
        if mod(t,1)==0
            posGPSNominal = trj.avp(k1,7:9)' + davp0(7:9).*randn(3,1);  % GPS with nominal white noise
            [posGPSFaulted, faultActive, ~] = applyGpsFaults(t, posGPSNominal, faultWindows, faultSpikes);
            yAll = insAll.pos - posGPSFaulted;
            yFD  = insFD.pos  - posGPSFaulted;
            kfAll = kfupdate(kfAll, yAll, 'M');
            kfFD  = kfupdate(kfFD,  yFD,  'M');
            [kfAll, insAll] = kffeedback(kfAll, insAll, 1, 'avp');
            [kfFD,  insFD]  = kffeedback(kfFD,  insFD,  1, 'avp');
            avpAll(ki,:)   = [insAll.avp', t];
            avpFD(ki,:)    = [insFD.avp',  t];
            xkpkAll(ki,:)  = [kfAll.xk; diag(kfAll.Pxk); t]';
            xkpkFD(ki,:)   = [kfFD.xk;  diag(kfFD.Pxk);  t]';
            nisLog(ki)         = kfFD.fd.nis;
            faultFlagLog(ki)   = faultActive;
            detectLog(ki)      = any(kfFD.fd.isOutlier);
            dwellLog(ki)       = any(kfFD.fd.inDwell);
            ki = ki + 1;
        end
        timebar;
    end

    avpAll(ki:end,:) = [];  avpFD(ki:end,:) = [];
    xkpkAll(ki:end,:) = []; xkpkFD(ki:end,:) = [];
    nisLog(ki:end) = []; faultFlagLog(ki:end) = [];
    detectLog(ki:end) = []; dwellLog(ki:end) = [];

    % Position error analysis
    tAll = avpAll(:,end);
    tFD  = avpFD(:,end);
    tMeas = tFD;                                  % time tags for logged measurement metrics
    truthAll = avpinterp1(trj.avp, tAll);
    truthFD  = avpinterp1(trj.avp, tFD);
    [RMhAll, clRNhAll] = RMRN(truthAll(:,7:9));
    [RMhFDTrue, clRNhFDTrue] = RMRN(truthFD(:,7:9));
    posErrAll = [(avpAll(:,7)-truthAll(:,7)).*RMhAll, ...
                 (avpAll(:,8)-truthAll(:,8)).*clRNhAll, ...
                 (avpAll(:,9)-truthAll(:,9))];
    posErrFD  = [(avpFD(:,7)-truthFD(:,7)).*RMhFDTrue, ...
                 (avpFD(:,8)-truthFD(:,8)).*clRNhFDTrue, ...
                 (avpFD(:,9)-truthFD(:,9))];
    horErrAll = sqrt(posErrAll(:,1).^2 + posErrAll(:,2).^2);
    horErrFD  = sqrt(posErrFD(:,1).^2 + posErrFD(:,2).^2);
    rmseAll  = sqrt(mean(posErrAll.^2,1));
    rmseFD   = sqrt(mean(posErrFD.^2,1));
    rmseAll3D = sqrt(mean(sum(posErrAll.^2,2)));
    rmseFD3D  = sqrt(mean(sum(posErrFD.^2,2)));

    % ANP computation from filter covariance (95% confidence)
    posVarFD = xkpkFD(:,15+(7:9));             % covariance of [dlat, dlon, dhgt]
    sigmaLat = sqrt(posVarFD(:,1));            % latitude 1-sigma in rad
    sigmaLon = sqrt(posVarFD(:,2));            % longitude 1-sigma in rad
    [RMhFDEst, clRNhFDEst] = RMRN(avpFD(:,7:9));   % meridian & transverse radii (m)
    sigmaNorth = RMhFDEst .* sigmaLat;         % convert to metres
    sigmaEast  = clRNhFDEst .* sigmaLon;
    kh95 = sqrt(5.9915);                       % sqrt(chi2inv(0.95,2)) for 95% circle
    anp = kh95 .* sqrt(sigmaNorth.^2 + sigmaEast.^2);

    sim = struct( ...
        'avpAll', avpAll, ...
        'avpFD', avpFD, ...
        'xkpkAll', xkpkAll, ...
        'xkpkFD', xkpkFD, ...
        'nisLog', nisLog, ...
        'faultFlagLog', faultFlagLog, ...
        'detectLog', detectLog, ...
        'dwellLog', dwellLog, ...
        'tAll', tAll, ...
        'tFD', tFD, ...
        'tMeas', tMeas, ...
        'posErrAll', posErrAll, ...
        'posErrFD', posErrFD, ...
        'horErrAll', horErrAll, ...
        'horErrFD', horErrFD, ...
        'rmseAll', rmseAll, ...
        'rmseFD', rmseFD, ...
        'rmseAll3D', rmseAll3D, ...
        'rmseFD3D', rmseFD3D, ...
        'anp', anp, ...
        'sigmaNorth', sigmaNorth, ...
        'sigmaEast', sigmaEast, ...
        'nisThreshold', kfFD.fd.chi2Threshold);
end


function [faultedPos, isFault, biasNEU] = applyGpsFaults(t, nominalPos, faultWindows, faultSpikes)
% Inject deterministic faults (in metres) into simulated GPS measurements.
    biasNEU = zeros(3,1);
    isFault = false;
    for k = 1:size(faultWindows,1)
        if t>=faultWindows(k,1) && t<faultWindows(k,2)
            biasNEU = biasNEU + faultWindows(k,3:5)';
            isFault = true;
        end
    end
    for k = 1:size(faultSpikes,1)
        if abs(t-faultSpikes(k,1))<1e-3
            biasNEU = biasNEU + faultSpikes(k,2:4)';
            isFault = true;
        end
    end
    if any(biasNEU)
        [RMh, clRNh] = RMRN(nominalPos');
        dPos = [biasNEU(1)/RMh(1); biasNEU(2)/clRNh(1); biasNEU(3)];
    else
        dPos = zeros(3,1);
    end
    faultedPos = nominalPos + dPos;
end

