function Tx_baseband = pA_to_basebandports(theta_steer_angle)
% This function is used to transmit the generated waveform 
% to the baseband ports of the transmitter (Tx)
% Input: theta angle / azimuth angle / steering direction
% Output: Tx_baseband signal: a 45x11,200 matrix

N = 45;             % Number of antennas
Tc = 28e-9;         % Clock frequecy
c = physconst('LightSpeed');
lambda = c*(1/Tc);  % wavelength
d = lambda/2;       % Inner-antenna spacing

%% Initial waveform generation (at Point A)

% Defined as a 1x11200 matrix

Amp = 1000;         % Amplitude is equal to 1kV
sig(1:7) = [-Amp -Amp -Amp Amp Amp -Amp Amp];
sig(8:1400) = zeros;  
signal = repmat(sig,1,8);

%% Use of weights
% Defined as a 45x1 matrix

w = ones(N,1);

Tx_baseband = w*signal; % 45x1 * 1x11,200 = 45x11,200 matrix

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

%% The original signal is phaseshifted

% Element-wise multiplication
Tx_baseband = v.*Tx_baseband;
   
end
