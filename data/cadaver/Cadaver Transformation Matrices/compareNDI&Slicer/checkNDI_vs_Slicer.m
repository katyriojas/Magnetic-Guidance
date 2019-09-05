% Katy Riojas
% Check readout from NDI track and Slicer saved files
% 8/22/2019
clear all; clc; close all;

ait_ndi_t = [9.12;-49.98;-1432.47];
ait_ndi_quat = [0.4526;0.3328;0.0559;0.8254]';
ait_ndi_R = quat2rotm(ait_ndi_quat);
T_ait_tracker_ndi = [ait_ndi_R,ait_ndi_t;0,0,0,1];

magsensor_ndi_t = [-249;230.13;-1177.64];
magsensor_ndi_quat = [0.4623;-0.4426;0.3116;-0.7023]';
magsensor_ndi_R = quat2rotm(magsensor_ndi_quat);
T_magsensor_tracker_ndi = [magsensor_ndi_R,magsensor_ndi_t;0,0,0,1];

mag_ndi_t = [27.86;43.65;-1476.70];
mag_ndi_quat = [0.1995;0.5170;0.1081;0.8253]';
mag_ndi_R = quat2rotm(mag_ndi_quat);
T_mag_tracker_ndi = [mag_ndi_R,mag_ndi_t;0,0,0,1];

fixture_ndi_t = [-21.91;-54.61;-1440.96];
fixture_ndi_quat = [0.2080;0.5302;0.1110;0.8145]';
fixture_ndi_R = quat2rotm(fixture_ndi_quat);
T_fixture_tracker_ndi = [fixture_ndi_R,fixture_ndi_t;0,0,0,1];

% Pull in Transforms
load('T_cochlea_tracker.mat');
T_fixture_tracker_slicer = AffineTransform_double_3_3;
load('T_ait_tracker.mat');
T_ait_tracker_slicer = AffineTransform_double_3_3;
load('T_mag_tracker.mat');
T_mag_tracker_slicer = AffineTransform_double_3_3;
load('T_magsensor_tracker.mat');
T_magsensor_tracker_slicer = AffineTransform_double_3_3;

% Convert
% Convert Transform to Tracker Space
T_fixture_tracker_slicer = savedTransform2TrackerSpace(T_fixture_tracker_slicer);
T_ait_tracker_slicer = savedTransform2TrackerSpace(T_ait_tracker_slicer);
T_mag_tracker_slicer = savedTransform2TrackerSpace(T_mag_tracker_slicer);
T_magsensor_tracker_slicer = savedTransform2TrackerSpace(T_magsensor_tracker_slicer);

% Differences
ait_diff = T_ait_tracker_ndi - T_ait_tracker_slicer
mag_diff = T_mag_tracker_ndi - T_mag_tracker_slicer
fixture_diff =  T_fixture_tracker_ndi - T_fixture_tracker_slicer
magsensor_diff =  T_magsensor_tracker_ndi - T_magsensor_tracker_slicer
