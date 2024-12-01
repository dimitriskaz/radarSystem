 function z_out = basebandports_to_pZ(Rx_baseband,theta_steer_angle)
% This function is used to pass the received signal 
% of the baseband ports Rx
% to the point Z
% Inputs:
%   1. Rx_baseband: the output of the backscatter function
%   2. theta angle / azimuth angle / steering direction
% Output:
%   1. signal at point Z (a 1x11,200 matrix)

N = 45;             % Number of antennas
Tc = 28e-9;         % Clock frequecy
c = physconst('LightSpeed');
lambda = c*(1/Tc);  % wavelength
d = lambda/2;       % Inner-antenna spacing

%% Use of phase-shifters

% Define r
% Use a more simplified code 
% (compared to the backscatterdata function)

r = (-N/2+1/2:1:N/2)*d;
r = cat(1,r,zeros(2,N));

% Derivation of u_theta and psi are explained in the report
u_theta = [cos(deg2rad(theta_steer_angle)),sin(deg2rad(theta_steer_angle)),0]';
psi = r'*(2*pi/lambda)*u_theta;
v = exp(-1i*psi);

%% The Rx baseband signal is phaseshifted

% Element-wise multiplication
Rx_baseband = conj(v).*Rx_baseband;

%% Use of weights

w = ones(N,1);
z_out = w'*Rx_baseband; % 1x45 * 45x11,200 = 1x11,200 matrix
 
end