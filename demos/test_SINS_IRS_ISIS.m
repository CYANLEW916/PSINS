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
