# Precise Strapdown Inertial Navigation System (PSINS) Toolbox for MATLAB

## 1. Preface

Precise Strapdown Inertial Navigation System (PSINS) toolbox for MATLAB is an open source program package, primarily developed for inertial-grade or higher grade inertial navigation system simulation and data processing. PSINS toolbox includes strapdown inertial sensor (gyro & accelerometer) sampling simulation, initial self-alignment simulation, pure SINS navigation algorithm simulation, SINS/DR & SINS/GPS integrated navigation simulation and many other useful routes, which are all implemented by a bunch of powerful library functions. The PSINS library functions are well modularized and organized, then they are easy to understand and master. Surely, PSINS toolbox has the capability to processing real SIMU and GPS sampling data with a little or even no modification. On the basis of this PSINS toolbox, users can quickly and conveniently set up an inertial navigation solution to achieve their specific purpose.

## 2. Symbols & Conventions

1) SINS: Strapdown Inertial Navigation System.
2) GPS: Global Positioning System.
3) DR: Dead Reckoning.
4) SIMU: Strapdown Inertial Measurement Unit.
5) POS: Positioning Orientation System.
6) i-frame: inertial frame.
7) e-frame: ECEF(Earth-Center-Earth-Fix) frame.
8) n-frame: navigation reference frame, with E-N-U (East-North-Up) pointing orientations.
9) b-frame: carrier’s body frame (i.e. SIMU frame), with X-Right, Y-Forward and Z-Up pointing orientations.
10) p-frame: computed navigation frame, misalignment angles from n-frame to p-frame usually denoted as ‘phi’.
11) ts/fs: SIMU sampling interval / sampling frequency.
12) T or len: total simulation time or data length.
13) att/qnb/Cnb: att – Euler’s angles of body attitude, i.e. att = [pitch, roll, yaw], NOTE: yaw angle follows the right-handed system convention with range (-pi, pi]; qnb – attitude quaternion representation; Cnb – direct cosine matrix (DCM), i.e. the transformation matrix from n-frame to b-frame.
14) vn: body velocity, i.e. the n-frame linear velocity w.r.t e-frame.
15) pos: body geographic coordinates, pos = [lat, lon, hgt], where lat – latitude, lon – longitude, hgt – height above sea level.
16) avp/avp0: avp = [att, vn, pos, t], t – time tag. Usually, avp0 specifies the initial navigation parameters, i.e. avp0 = [att0, vn0, pos0].
17) eb/web/db/wdb: eb – gyro constant drift error; web – gyro angular random walk coefficient; db – accelerometer constant bias; wdb – accelerometer velocity random walk coefficient.
18) taug/taua: correlation time for gyro/accelerometer 1st order Markov process.
19) dKg/dKa: scale factor errors and misalignment errors for gyro/accelerometer triad.
20) imuerr: structure array including eb, web, db, wdb, taug, taua, dKg and dKa.
21) wm/vm: the increment of gyro angular/accelerometer velocity sampling data within ts. Sometimes, symbols wib/fb are used for gyro angular rate / accelerometer specific force.
22) imu: imu = [wm, vm, t], t – time tag, each row of imu represents a SIMU incremental sample within [t-ts,t].
23) wnie/wnen/wnin: wnie – the Earth’s angular rate projection in n-frame; wnen – the rotation rate due to body’s linear motion on the Earth’s surface; wnin = wnie+wnen.
24) gn: gravity vector.
25) phi/dvn/dpos: phi – platform errors (small misalignment angles) from n-frame to p-frame; dvn – velocity errors; dpos – geographic position errors, dpos = [dlat, dlon, dhgt].
26) davp/davp0: davp = [phi, dvn, dpos, t]. Usually, davp0 specifies the initial avp error, i.e. davp0 = [phi0, dvn0, dpos0].
27) eth: structure array, containing some important/useful parameters related to the Earth’s navigation model.
28) trj: trajectory simulation result, including SIMU sensor outputs and trajectory avp references, etc.
29) ins: structure array for SINS updating algorithm.
30) kf: Kalman filter structure array.
31) xkpk: the results of Kalman filter updating, including state estimation, diagonal of covariance matrix and time tag, i.e. xkpk = [kf.xk, diag(kf.Pxk), t].

## 3. System Requirements

When developing this toolbox, the author’s PC setting is: 

Microsoft Windows 7 (SP1) + MATLAB 8.2.0 (R2013b) + CPU 2.1GHz + RAM 2.0GB.

## 4. Quick Start

1. Copy the PSINS toolbox root folder `psins\`, including all subfolders and files, to your computer.
2. Run `psins\psinsinit.m` to initialize PSINS environment.
3. Run `psins\demos\test_SINS_trj.m` to generate a moving trajectory.
4. Run `psins\demos\test_SINS_GPS.m` to demonstrate SINS/GPS integrated navigation.
5. There are many demo examples in `psins\demos`, such as coning & sculling motion demonstration, initial alignment, pure inertial navigation and POS data fusion, etc.
6. Try to do some modification and put your exercise file under `psins\mytest`. Enjoy yourself and may you find something helpful!

