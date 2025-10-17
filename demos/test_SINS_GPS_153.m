% SINS/GPS intergrated navigation simulation using kalman filter.
% Please run 'test_SINS_trj.m' to generate 'trj10ms.mat' beforehand!!!
% See also  test_SINS_trj, test_SINS, test_SINS_GPS_186, test_SINS_GPS_193.
% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 17/06/2011
glvs
psinstypedef(153);
trj = trjfile('trj10ms.mat');
% initial settings
[nn, ts, nts] = nnts(2, trj.ts);
imuerr = imuerrset(0.03, 100, 0.001, 5);
imu = imuadderr(trj.imu, imuerr);
davp0 = avperrset([0.5;-0.5;20], 0.1, [1;1;3]);
ins = insinit(avpadderr(trj.avp0,davp0), ts);
% KF filter
rk = poserrset([1;1;3]);
kf = kfinit(ins, davp0, imuerr, rk);
kf.Pmin = [avperrset(0.01,1e-4,0.1); gabias(1e-3, [1,10])].^2;  kf.pconstrain=1;
len = length(imu); [avp, xkpk] = prealloc(fix(len/nn), 10, 2*kf.n+1);
timebar(nn, len, '15-state SINS/GPS Simulation.'); 
ki = 1;
for k=1:nn:len-nn+1
    k1 = k+nn-1;  
    wvm = imu(k:k1,1:6);  t = imu(k1,end);
    ins = insupdate(ins, wvm);
    kf.Phikk_1 = kffk(ins);
    kf = kfupdate(kf);
    if mod(t,1)==0
        posGPS = trj.avp(k1,7:9)' + davp0(7:9).*randn(3,1);  % GPS pos simulation with some white noise
        kf = kfupdate(kf, ins.pos-posGPS, 'M');
        [kf, ins] = kffeedback(kf, ins, 1, 'avp');
        avp(ki,:) = [ins.avp', t];
        xkpk(ki,:) = [kf.xk; diag(kf.Pxk); t]';  ki = ki+1;
    end
    timebar;
end
avp(ki:end,:) = [];  xkpk(ki:end,:) = []; 
% show results
insplot(avp);
avperr = avpcmpplot(trj.avp, avp);
kfplot(xkpk, avperr, imuerr);

% Actual Navigation Performance (ANP) estimation and RNP comparison
global glv
posVar = xkpk(:,15+(7:9));                 % covariance of [dlat, dlon, dhgt]
sigmaLat = sqrt(posVar(:,1));               % latitude 1-sigma in rad
sigmaLon = sqrt(posVar(:,2));               % longitude 1-sigma in rad
[RMh, clRNh] = RMRN(avp(:,7:9));            % meridian & transverse radii (m)
sigmaNorth = RMh .* sigmaLat;               % convert to metres
sigmaEast = clRNh .* sigmaLon;
kh95 = sqrt(5.9915);                        % sqrt(chi2inv(0.95,2)) for 95% circle
anp = kh95 .* sqrt(sigmaNorth.^2 + sigmaEast.^2);
rnp = 0.1 * glv.nm * ones(size(anp));       % RNP AR 0.1 requirement in metres

myfigure('ANP_vs_RNP');
plot(avp(:,end), anp, 'b', avp(:,end), rnp, 'r--', 'LineWidth', 1.5);
hold on;
violations = anp > rnp;
if any(violations)
    plot(avp(violations,end), anp(violations), 'ro', 'MarkerFaceColor', 'r');
end
grid on;  xygo('Time / s', 'Horizontal performance / m');
legend('ANP (95%)','RNP AR 0.1','Location','best');


