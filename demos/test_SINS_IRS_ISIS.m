% SINS sensor data simulation for IRS1/2 and ISIS based on test_SINS.m.
% Requires 'trj_NLG_approach.mat' from NLG approach trajectory data.
% See also  test_SINS, test_SINS_trj, imuerrset, imuadderr, inspure.

glvs
trjPath = 'D:\All Model\psins251010\data\trj_NLG_approach.mat';
trj = trjfile(trjPath);

%% initial AVP error setting
initAvpErr = avperrset([0.5;0.5;5], 0.1, [10;10;10]);
avp00 = avpadderr(trj.avp0, initAvpErr);
trj = bhsimu(trj, 1, 10, 3, trj.ts);

%% sensor specs and simulation
specs = makeSensorSpecs();
sensorData = repmat(struct('name','', 'imu',[], 'avp',[]), numel(specs), 1);

for k = 1:numel(specs)
    rng(100 + k, 'twister');
    imuerr = imuerrset(specs(k).gyroBias, specs(k).accBias, ...
        specs(k).gyroArw, 0);
    imu = imuadderr(trj.imu, imuerr);
    avp = inspure(imu, avp00, trj.bh, 1);
    avp = addOutputNoise(avp, specs(k).attNoise, specs(k).posNoise);

    sensorData(k).name = specs(k).name;
    sensorData(k).imu = imu;
    sensorData(k).avp = avp;
end

save('trj10ms_sensor_data.mat', 'sensorData', 'specs', 'trj');

%% helper functions
function specs = makeSensorSpecs()
%MAKESENSORSPECS Build sensor spec table for IRS1/2 and ISIS.
%   specs = MAKESENSORSPECS() returns a struct array with error specs
%   (gyro bias/ARW, accel bias) and output noise (attitude/position).
    specs = struct('name', {'INS1', 'INS2', 'ISIS'}, ...
        'gyroBias', {0.01, 0.01, 0.1}, ...
        'gyroArw',  {0.003, 0.003, 0.1}, ...
        'accBias',  {100, 100, 5000}, ...
        'attNoise', {0, 0, 0}, ...
        'posNoise', {0, 0, 0});
end

function avp = addOutputNoise(avp, attNoiseDeg, posNoiseM)
%ADDOUTPUTNOISE Add attitude/position output noise to AVP.
%   avp = ADDOUTPUTNOISE(avp, attNoiseDeg, posNoiseM) applies 1-sigma
%   attitude noise (deg) and position noise (m) to the AVP sequence.
    global glv
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
