%% Task 2 

function Rx_baseband = backscatterdata(Tx_baseband, targets)
% This function is used to construct the channel

% Inputs:
%   1. Tx_baseband is defined outside the function (it is a 45xsnapshots matrix
%   determined by the pA_to_basebandports function)
%   2. Parameter target depends on the number of targets set for each Task
%   Each target is defined as by the parameters: theta, R, RCS

% Outputs:
%   1. Rx_baseband signal(it is matrix of 45xsnapshots = 45x11,200)

N = 45;             % Number of array of antennas
Tc = 28e-9;         % Clock fr
fc = 15e9;          % Carrier frequency at the middle of Ku-band
c = physconst('LightSpeed');
lambda = c/(fc);    % Wavelength
d = lambda/2;       % Inner-antenna spacing
Tp = 7*Tc;          % Pulse duration
Bw = 1/Tc;          % Bandwidth
one_PRI = 200*Tp;   % Pulse repetition
numPRI = 8;         % Number of PRIs
%Pfa = 0.1;         % Probability of false alarm
Ptx = 1;            % Power of the transmitter Tx
Gtx = 1;
Grx = 1;
M = 199;            % Bins
snapshots = (M+1)*7*8;
kB = 1.28e-23;      % Boltzman constant
Temp = 290;         % Temperature
Fn = 3.1068;        % Noise factor

% Initialisatiion of Rx_baseband: the output of the function
Rx_baseband = zeros(N,snapshots);

% Define r
r = [-22*d,-21*d,-20*d,-19*d,-18*d,-17*d,-16*d,-15*d,-14*d,-13*d,-12*d,-11*d,-10*d,-9*d,-8*d,-7*d,-6*d,-5*d,-4*d,-3*d,-2*d,-d,0,d,2*d,3*d,4*d,5*d,6*d,7*d,8*d,9*d,10*d,11*d,12*d,13*d,14*d,15*d,16*d,17*d,18*d,19*d,20*d,21*d,22*d;...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

for j=1:1:targets
    
    %% Define the target's parameters
    
    RCS = [];
    ns = 1;
    
    for PRI=1:1:numPRI
        if targets == 1
            theta_steer_angle = 40;
            R = 2000;
            RCS_mean = 1;
            % Constant RCS
            RCS = [RCS RCS_mean*ones(ns,snapshots/numPRI)];
        end
        if targets == 2
            theta_steer_angle = 70;
            R = 3000;
            RCS_mean = 5;
            % Swerlingmodel1&2 model: scatters of similar amplitudes
            RCS = [RCS exprnd(RCS_mean, [ns,ns])*ones(ns,snapshots/numPRI)];
        end
        if targets == 3
            theta_steer_angle = 120;
            R = 2500;
            RCS_mean = 4.5;
            % Swerlingmodel3&4 model: scatters with one much larger than
            % the other
            RCS = [RCS (RCS_mean/4)*chi2rnd(4, [ns,ns])*ones(ns,snapshots/numPRI)];
        end
    end
    
    %% Tx array manifold vectors
    
    phi_steer_angle = 0; % elevation angle
    u_theta = [cos(deg2rad(theta_steer_angle))*cos(deg2rad(phi_steer_angle));...
        sin(deg2rad(theta_steer_angle))*cos(deg2rad(phi_steer_angle));...
        sin(deg2rad(phi_steer_angle))];
    Stx = exp(1i*2*pi*r'*u_theta/lambda); % 45 x 1 vector
    Tx_out = Stx.'*Tx_baseband; % 1x45 * 45x11,200 = 1x11,200 vector
    
     %% Define delay
    
    %delay = 0; % for testing (see Appendix A)
    techo = 2*R/c;
    delay = techo/Tc;
    % rounds the elements to the nearest integer 
    % greater than or equal to that element
    delay = ceil(delay);
    
    %% Path Attenuation
    
    %betas = 1; % for testing (see Appendix A)
    betas = sqrt(Ptx*Gtx*Grx/(4*pi)^3)*(lambda/R^2)*sqrt(RCS);
    Tx_out2 = betas.*Tx_out; % 1x11,200 vector
    
    %% Received signal
    
    Srx = exp(-1i*2*pi*r'*u_theta/lambda); % 45 x 1 vector
    Rx_out = Srx*Tx_out2; % 45x1 * 1x11,200 = 45x11,200 vector
    
    %% Backscatterdata at basenand ports of Rx
    
    Rx_baseband(:,1+delay:snapshots) = Rx_baseband(:,1+delay:snapshots) + Rx_out(:,1:end-delay);
    
end

%% Noise

% Noise: normally distributed random numbers
noise = sqrt(kB*Temp*Fn*Bw)*(randn(N,snapshots)+1i*randn(N,snapshots))/sqrt(2); 
Rx_baseband = Rx_baseband + noise;

end