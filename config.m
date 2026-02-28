function cfg = config()
% CONFIG  Parameter configuration for PPV-aided SWGLT fault detection.
%   cfg = CONFIG() returns a structure containing all simulation parameters
%   including sensor noise, fault injection, algorithm settings, and plotting.
%
% See also  main_simulation, fdi_glt, fdi_wglt, fdi_swglt.

    %% Simulation parameters
    cfg.ts = 0.02;           % sampling period (s), 50 Hz
    cfg.rng_seed = 42;       % random seed for reproducibility

    %% Sensor noise characteristics (from paper Table 1)
    % Unit conversion constants
    dph2rps = pi / 180 / 3600;  % 1 deg/h in rad/s
    g0 = 9.8;                   % gravitational acceleration (m/s^2)

    % Gyroscope bias instability (deg/h -> rad/s)  [per DO-236E nav-grade IRS]
    cfg.sigma_ins_gyro  = 0.005 * dph2rps;  % INS1 & INS2 gyro noise std (rad/s)
    cfg.sigma_isis_gyro = 0.05  * dph2rps;  % ISIS gyro noise std (rad/s)

    % Accelerometer bias instability (g -> m/s^2)  [per DO-236E nav-grade IRS]
    cfg.sigma_ins_acc  = 5e-5 * g0;         % INS1 & INS2 accel noise std (m/s^2)
    cfg.sigma_isis_acc = 1e-3 * g0;         % ISIS accel noise std (m/s^2)

    % Noise standard deviation vectors (9 sensors: INS1 xyz, INS2 xyz, ISIS xyz)
    cfg.sigma_gyro = [cfg.sigma_ins_gyro  * ones(1,3), ...
                      cfg.sigma_ins_gyro  * ones(1,3), ...
                      cfg.sigma_isis_gyro * ones(1,3)];  % 1x9

    cfg.sigma_acc  = [cfg.sigma_ins_acc  * ones(1,3), ...
                      cfg.sigma_ins_acc  * ones(1,3), ...
                      cfg.sigma_isis_acc * ones(1,3)];   % 1x9

    %% Configuration matrix H (9x3)
    % H = kron(ones(3,1), eye(3)) or equivalently repmat(eye(3), 3, 1)
    cfg.H = repmat(eye(3), 3, 1);
    cfg.n_sensors = 9;      % total number of sensors per type
    cfg.n_dof = 3;          % degrees of freedom (3 axes)
    cfg.n_parity = cfg.n_sensors - cfg.n_dof;  % parity space dimension = 6

    %% Fault injection parameters
    % Condition 1: INS2 Z-axis gyro hard fault (step fault)
    cfg.fault1.sensor_idx = 6;          % INS2 Z-axis (6th element in Z_gyro)
    cfg.fault1.type = 'gyro';
    cfg.fault1.intervals = [120 135; 440 455; 595 610];  % time intervals (s)
    cfg.fault1.magnitudes = [0.8; 1.5; 3.0] * dph2rps;  % fault magnitudes (rad/s)

    % Condition 2: INS1 Y-axis accelerometer hard fault
    cfg.fault2.sensor_idx = 2;          % INS1 Y-axis (2nd element in Z_acc)
    cfg.fault2.type = 'acc';
    cfg.fault2.intervals = [100 115; 420 435; 575 590];  % time intervals (s)
    cfg.fault2.magnitudes = [0.008; 0.015; 0.03] * g0;  % fault magnitudes (m/s^2)

    % Condition 3: INS1 X-axis gyro soft fault (ramp fault)
    cfg.fault3.sensor_idx = 1;          % INS1 X-axis (1st element in Z_gyro)
    cfg.fault3.type = 'gyro';
    cfg.fault3.interval = [260 310];    % time interval (s)
    cfg.fault3.rate = 0.015 * dph2rps;  % fault rate 0.015 deg/h per second
    cfg.fault3.t_start_ref = 260;       % reference start time for ramp: fault = rate*(t-260)

    %% GLT parameters
    cfg.PFA = 0.01;                      % probability of false alarm
    cfg.T_D = chi2inv(1 - cfg.PFA, cfg.n_parity);  % detection threshold ~16.81

    %% SWGLT parameters
    cfg.window_length = 25;              % PCA window = 0.5s / 0.02s = 25 samples
    cfg.pca_threshold = 0.85;            % cumulative variance ratio for PCA
    cfg.K_tolerance = 0.01;              % tolerance factor for adaptive threshold

end
