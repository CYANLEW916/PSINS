% Fault detection demo for IRS1/2 and ISIS sensor outputs.
% Requires 'trj10ms_sensor_data.mat' from test_SINS_IRS_ISIS.m.
% See also  test_SINS_IRS_ISIS, avpcmpplot, chi2inv.

glvs
load('trj10ms_sensor_data.mat', 'sensorData', 'specs', 'trj');

%% residuals (AVP)
idxIrs1 = find(strcmp({sensorData.name}, 'IRS1'), 1);
idxIrs2 = find(strcmp({sensorData.name}, 'IRS2'), 1);
idxIsis = find(strcmp({sensorData.name}, 'ISIS'), 1);

avp1 = sensorData(idxIrs1).avp;
avp2 = sensorData(idxIrs2).avp;
avpI = sensorData(idxIsis).avp;

% Align lengths with reference trajectory (bhsimu may change length).
n = min([size(avp1, 1), size(avp2, 1), size(avpI, 1), size(trj.avp, 1)]);
avp1 = avp1(1:n, :);
avp2 = avp2(1:n, :);
avpI = avpI(1:n, :);
trjAvp = trj.avp(1:n, :);

r12 = avp1(:, 1:9) - avp2(:, 1:9);
r1I = avp1(:, 1:9) - avpI(:, 1:9);
r2I = avp2(:, 1:9) - avpI(:, 1:9);

r1T = avp1(:, 1:9) - trjAvp(:, 1:9);
r2T = avp2(:, 1:9) - trjAvp(:, 1:9);

%% chi-square thresholds (attitude + position)
alpha = 0.997; % about 3-sigma
[R12, R1I] = buildResidualCov(specs(idxIrs1), specs(idxIsis));
thr12 = chi2inv(alpha, size(R12, 1));
thr1I = chi2inv(alpha, size(R1I, 1));

j12 = rowChiSquare(r12, R12);
j1I = rowChiSquare(r1I, R1I);
j2I = rowChiSquare(r2I, R1I);

flag12 = j12 > thr12;
flag1I = j1I > thr1I;
flag2I = j2I > thr1I;

%% simple isolation logic
faultLabel = strings(size(j12));
for k = 1:numel(j12)
    if flag12(k)
        if flag1I(k) && flag2I(k)
            faultLabel(k) = isolateByTruth(r1T(k, :), r2T(k, :));
        else
            faultLabel(k) = "IRS1_or_IRS2";
        end
    elseif flag1I(k) && flag2I(k)
        faultLabel(k) = "ISIS";
    else
        faultLabel(k) = "Normal";
    end
end

%% plots
figure; subplot(3, 1, 1);
plot(j12); hold on; yline(thr12, 'r--'); grid on;
title('IRS1-IRS2 Chi-square');
subplot(3, 1, 2);
plot(j1I); hold on; yline(thr1I, 'r--'); grid on;
title('IRS1-ISIS Chi-square');
subplot(3, 1, 3);
plot(j2I); hold on; yline(thr1I, 'r--'); grid on;
title('IRS2-ISIS Chi-square');

%% helper functions
function [R12, R1I] = buildResidualCov(irsSpec, isisSpec)
%BUILDRESIDUALCOV Build residual covariance for attitude/position.
%   [R12, R1I] = BUILDRESIDUALCOV(irsSpec, isisSpec) returns covariance
%   matrices for IRS1-IRS2 and IRS-ISIS residuals.
    global glv
    attIrs = irsSpec.attNoise * glv.deg;
    posIrs = irsSpec.posNoise;
    attIsis = isisSpec.attNoise * glv.deg;
    posIsis = isisSpec.posNoise;

    sigma12 = [attIrs * sqrt(2) * ones(3, 1); zeros(3, 1); ...
        posIrs * sqrt(2) * ones(3, 1)];
    sigma1I = [sqrt(attIrs^2 + attIsis^2) * ones(3, 1); zeros(3, 1); ...
        sqrt(posIrs^2 + posIsis^2) * ones(3, 1)];

    sigma12(7:8) = sigma12(7:8) / glv.Re;
    sigma1I(7:8) = sigma1I(7:8) / glv.Re;

    R12 = diag(sigma12 .^ 2);
    R1I = diag(sigma1I .^ 2);
end

function jval = rowChiSquare(residual, R)
%ROWCHISQUARE Row-wise chi-square statistic.
%   jval = ROWCHISQUARE(residual, R) returns J for each row.
    rinv = inv(R);
    jval = sum((residual * rinv) .* residual, 2);
end

function label = isolateByTruth(r1T, r2T)
%ISOLATEBYTRUTH Decide which IRS deviates more from truth.
%   label = ISOLATEBYTRUTH(r1T, r2T) returns "IRS1" or "IRS2".
    if norm(r1T) >= norm(r2T)
        label = "IRS1";
    else
        label = "IRS2";
    end
end
