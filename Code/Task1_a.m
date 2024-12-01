%% Task 1: (a) Phase Shifters
% The following code estimates the vectors of Tx and Rx phase-shift
% for any steering direction theta (defined as theta_steer_angle)
clc;
clear;
close all

N = 45;                             % Number of array of antennas
Tc = 28e-9;                         % clock frequency
c = physconst('LightSpeed');
lambda = c*Tc;                      % Wavelength
d = lambda/2;

% Defining r
r = [-22*d,-21*d,-20*d,-19*d,-18*d,-17*d,-16*d,-15*d,-14*d,-13*d,-12*d,-11*d,-10*d,-9*d,-8*d,-7*d,-6*d,-5*d,-4*d,-3*d,-2*d,-d,0,d,2*d,3*d,4*d,5*d,6*d,7*d,8*d,9*d,10*d,11*d,12*d,13*d,14*d,15*d,16*d,17*d,18*d,19*d,20*d,21*d,22*d;...
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
     
theta = 0:1:360;                    % azimuth
phi = -90:1:90;                     % elevation   
% Calculation of vector by setting the direction in degrees
theta_steer_angle = 40;            % azimuth angle
phi_steer_angle = 0;                % elevation angle
k_theta = [cos(deg2rad(theta_steer_angle))*cos(deg2rad(phi_steer_angle));...
          sin(deg2rad(theta_steer_angle))*cos(deg2rad(phi_steer_angle));...
          sin(deg2rad(phi_steer_angle))];
w = exp(-1i*2*pi/lambda*r'*k_theta);
phase_shifters = angle(conj(w));
phase_shifters = rad2deg(phase_shifters);
% The phase shifter vector
phase_shifters_deg = mod(phase_shifters,360);