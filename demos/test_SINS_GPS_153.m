% SINS/GPS integrated navigation simulation with residual-based fault detection.
% Please run 'test_SINS_trj.m' to generate 'trj10ms.mat' beforehand!!!
% See also  test_SINS_trj, test_SINS, test_SINS_GPS_186, test_SINS_GPS_193.
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 17/06/2011, updated 15/11/2024
glvs
psinstypedef(153);
trj = trjfile('trj10ms.mat');

% initial settings
[nn, ts, nts] = nnts(2, trj.ts);
imuerr = imuerrset(0.05, 200, 0.002, 10);
imu = imuadderr(trj.imu, imuerr);
davp0 = avperrset([0.5;-0.5;20], 0.1, [1;1;3]);
insAll = insinit(avpadderr(trj.avp0,davp0), ts);
insFD  = insinit(avpadderr(trj.avp0,davp0), ts);

% Kalman filters
gpsMeasStdENU = [4; 4; 7];               % GPS measurement 1-sigma in metres (N, E, U)
rk = poserrset(gpsMeasStdENU);
kfAll = kfinit(insAll, davp0, imuerr, rk);
kfFD  = kfinit(insFD,  davp0, imuerr, rk);
commonPmin = [avperrset(0.01,1e-4,0.1); gabias(1e-3, [1,10])].^2;
kfAll.Pmin = commonPmin;  kfAll.pconstrain = 1;
kfFD.Pmin  = commonPmin;  kfFD.pconstrain  = 1;

% residual-based fault detection configuration (NIS gate with 99% chi-square threshold)
kfFD.fd.enable = 1;
kfFD.fd.chi2Threshold = 11.345;  % ~chi2inv(0.99, 3)
kfFD.fd.holdTime = 5;            % dwell time (s) before re-enabling rejected GPS measurements

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

gpsNoiseStd = poserrset(gpsMeasStdENU);   % convert to latitude/longitude radians & height metres

len = length(imu);
rows = fix(len/nn);
[avpAll, avpFD, xkpkAll, xkpkFD, nisLog, faultFlagLog, detectLog, dwellLog] = ...
    prealloc(rows, 10, 10, 2*kfAll.n+1, 2*kfFD.n+1, 1, 1, 1, 1);
timebar(nn, len, '15-state SINS/GPS Simulation with residual FDI.');
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
        posGPSNominal = trj.avp(k1,7:9)' + gpsNoiseStd.*randn(3,1);  % GPS with increased white noise
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

% show results
insplot(avpFD);  % best-performing filter
avpcmpplot(trj.avp, avpAll, avpFD);
avperrFD = avpcmp(avpFD, trj.avp);
kfplot(xkpkFD, avperrFD, imuerr);

% Position error analysis
global glv
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
fprintf('\nPosition RMSE (All measurements) [North East Up] / m: %6.2f %6.2f %6.2f\n', rmseAll);
fprintf('Position RMSE (FD-enabled)      [North East Up] / m: %6.2f %6.2f %6.2f\n', rmseFD);
fprintf('3D position RMSE - All: %6.2f m, FD-enabled: %6.2f m\n\n', rmseAll3D, rmseFD3D);

myfigure('Position Error Comparison');
labels = {'North', 'East', 'Up'};
for idx = 1:3
    subplot(3,1,idx);
    plot(tAll, posErrAll(:,idx), 'b', tFD, posErrFD(:,idx), 'r--', 'LineWidth', 1.2);
    grid on;  ylabel([labels{idx}, ' / m']);
    if idx==1
        title('Position errors with and without residual-based fault rejection');
    end
    if idx==3
        xlabel('Time / s');
    end
    legend('All measurements','FD-enabled','Location','best');
end

% residual-based fault detection metrics
nisThreshold = kfFD.fd.chi2Threshold;
fdActive = detectLog | dwellLog;
myfigure('Residual-based Fault Detection');
subplot(2,1,1);
plot(tMeas, nisLog, 'b', 'LineWidth', 1.2); hold on;
yline(nisThreshold, 'r--', 'LineWidth', 1.2);
grid on;  xlabel('Time / s');  ylabel('NIS');
legend('NIS','Threshold','Location','best');
title('Normalized innovation squared (NIS) test');
subplot(2,1,2);
stairs(tMeas, faultFlagLog, 'k', 'LineWidth', 1.2); hold on;
stairs(tMeas, fdActive, 'r--', 'LineWidth', 1.2);
stairs(tMeas, dwellLog, 'Color',[0.3 0.6 0.9], 'LineWidth', 1.1);
grid on;  xlabel('Time / s');
ylabel('State');
ylim([-0.1 1.1]);
legend('Injected fault','FD decision','Dwell active','Location','best');

% Actual Navigation Performance (ANP) estimation versus realised horizontal error (FD-enabled)
posVarFD = xkpkFD(:,15+(7:9));             % covariance of [dlat, dlon, dhgt]
sigmaLat = sqrt(posVarFD(:,1));            % latitude 1-sigma in rad
sigmaLon = sqrt(posVarFD(:,2));            % longitude 1-sigma in rad
[RMhFDEst, clRNhFDEst] = RMRN(avpFD(:,7:9));   % meridian & transverse radii (m)
sigmaNorth = RMhFDEst .* sigmaLat;         % convert to metres
sigmaEast  = clRNhFDEst .* sigmaLon;
kh95 = sqrt(5.9915);                       % sqrt(chi2inv(0.95,2)) for 95% circle
anp = kh95 .* sqrt(sigmaNorth.^2 + sigmaEast.^2);

myfigure('ANP_vs_HorizontalError');
plot(tFD, anp, 'b', 'LineWidth', 1.5);
hold on;
plot(tFD, horErrFD, 'Color',[0.1 0.6 0.1], 'LineWidth', 1.5);
rnpAR01 = 92.6;                           % RNP AR 0.1 NM requirement (metres)
yline(rnpAR01, 'r-', 'LineWidth', 1.2);
errBeyondANP = horErrFD > anp;
if any(errBeyondANP)
    plot(tFD(errBeyondANP), horErrFD(errBeyondANP), 'ks', 'MarkerFaceColor', [0.1 0.6 0.1]);
end
grid on;  xygo('Time / s', 'Horizontal performance / m');
legend('ANP (95%)','Horizontal error','RNP AR 0.1 (92.6 m)','Location','best');


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

