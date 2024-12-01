%% Task 5: 3rd scan - two targets

clc;
clear;
close all

%% Part a: Target-1 and Target-2 parameters are known: generate synthetic backscatter data for this scan
for targets = 1:1:2
i = 1;
backscatterData = cell(121,1);
for theta_steer_angle = 30:1:150
    
    % Define signal at the baseband ports of Tx
    Tx_baseband = pA_to_basebandports(theta_steer_angle);
    
    % Generate the backscatter data at the baseband ports of Rx
    Rx_baseband = backscatterdata(Tx_baseband, targets);
   
    backscatterData{i,1} = Rx_baseband;
    i = i+1;
end

%% Part b: Plot the backscatter-data for the dwell time that corresponds to the direction of Target-1 and Target-2

z_out = [];

for i = 1:1:121
    Rx_baseband = cell2mat(backscatterData([i]));
    theta_steer_angle = i+29;
    z_out(i,:) = basebandports_to_pZ(Rx_baseband,theta_steer_angle);
end

if targets == 1
    theta_steer_angle1 = 40;
    j1 = theta_steer_angle1 - 29;
    z1 = z_out(j1,:);
end

if targets == 2
    theta_steer_angle2 = 70;
    j2 = theta_steer_angle2 - 29;
    z2 = z_out(j2,:);
end


end 

figure();
semilogy(abs(z1));
title('Signal at point Z due to Target-1');
xlabel('Time Index (Tc)');
ylabel('Magnitude (V)');
xlim([0 11200]);
ylim([1e-7 1e-3]);
hold on;
Nc = 7;
for i=1:1:(Nc+1)
    xline(1400*i, 'red', ['PRI = ', num2str(i)], 'LineWidth', 1,'LabelHorizontalAlignment', 'Left');
end
set(gca, 'Fontsize', 14);

figure();
semilogy(abs(z2));
title('Signal at point Z due to Target-2');
xlabel('Time Index (Tc)');
ylabel('Magnitude (V)');
xlim([0 11200]);
ylim([1e-7 1e-3]);
hold on;
Nc = 7;
for i=1:1:(Nc+1)
    xline(1400*i, 'red', ['PRI = ', num2str(i)], 'LineWidth', 1,'LabelHorizontalAlignment', 'Left');
end
set(gca, 'Fontsize', 14);

% figure();
% semilogy(abs(z1));
% hold on;
% semilogy(abs(z2));
% title('Signal at point Z due to Target-1 and Target-2');
% xlabel('Time Index (Tc)');
% ylabel('Magnitude (V)');
% xlim([0 11200]);
% ylim([1e-7 1e-3]);
% hold on;
% Nc = 7;
% for i=1:1:(Nc+1)
%     xline(1400*i, 'red', ['PRI = ', num2str(i)], 'LineWidth', 1,'LabelHorizontalAlignment', 'Left');
% end
% legend('Target-1' , 'Target-2', 'Location', 'Northeast');
% set(gca, 'Fontsize', 14);

%% Part c: Using only your generated backscatter random numbers detect and estimate the parameters of the two targets

%% Define parameters required

Tc = 28e-9;
Fc = 15e9;
c = 3e8;
lambda = c/Fc;
Nc = 7;
M = 199;
Ptx = 1;

for targets = 1:1:2
    
%% Find the noise power at point Z

% Set a random theta_angle
theta_steer_angle = 30;

% Define signal at the baseband ports of Tx
Tx_baseband = pA_to_basebandports(theta_steer_angle);

% Define signal at the baseband ports of Rx
% Zero targets are used
Rx_baseband = backscatterdata(Tx_baseband, 0);

% Generate noise signal at point Z
z_noise = basebandports_to_pZ(Rx_baseband,theta_steer_angle);
zn = abs(z_noise);

noise_sq = zn.^2; 
% Estimate the noise power
noise_power = mean(noise_sq);

%% Find the threshold

% Define the pdf of the noise
% Add the threshold
% The following code uses a part of the suggested code in Class 4

vol_max = 1e-5;
vol_step = 1e-6;
voltage = 0:vol_step:vol_max;
noise_scale = sqrt(noise_power/(4-pi)/2); 
noise = raylpdf(voltage, noise_scale);
% Probability of false alarm set at 10^-3
Pfa = 0.001; 
y = cumtrapz(noise);
x = cumtrapz(noise)*vol_step;
vol_index = find(cumtrapz(noise)*vol_step<(1-Pfa));
v_thres = voltage(max(vol_index));

max_val = v_thres;
threshold = v_thres;

% Define an initial value for the angle that is detected
target_angle = 0;

%% Scan and find the steering direction of the target

if targets == 1
    i = 1;
    backscatterData = cell(121,1);
    for theta_steer_angle = 30:1:150
        Tx_baseband = pA_to_basebandports(theta_steer_angle);
        Rx_baseband = backscatterdata(Tx_baseband, targets);
        backscatterData{i,1} = Rx_baseband;
        i = i+1;
    end
    
    z_out = [];
    
    for i = 1:1:121
        Rx_baseband = cell2mat(backscatterData([i]));
        theta_steer_angle = i+29;
        z_out(i,:) = basebandports_to_pZ(Rx_baseband,theta_steer_angle);
    end
    
    z = abs(z_out); 

    for i = 1:1:121
        maximum = max(z(i,:));
        if maximum > max_val
            max_val = maximum;
            target_angle = i;
        end
    end
end

if targets == 2
    i = 1;
    backscatterData = cell(121,1);
    for theta_steer_angle = 30:1:150
        Tx_baseband = pA_to_basebandports(theta_steer_angle);
        Rx_baseband = backscatterdata(Tx_baseband, targets);
        backscatterData{i,1} = Rx_baseband;
        i = i+1;
    end
    
    z_out = [];
    
    for i = 1:1:121
        Rx_baseband = cell2mat(backscatterData([i]));
        theta_steer_angle = i+29;
        z_out(i,:) = basebandports_to_pZ(Rx_baseband,theta_steer_angle);
    end
    
    z = abs(z_out); 

    for i = 1:1:121
        maximum = max(z(i,:));
        if maximum > max_val
            max_val = maximum;
            target_angle = i;
        end
    end
end


%% Setting the parameters for the target according to the steering direction found

target_found = target_angle + 29;

% Define signal at the baseband ports of Tx
Tx_baseband = pA_to_basebandports(target_found);

% Define signal at the baseband ports of Rx
Rx_baseband = backscatterdata(Tx_baseband,targets);

% Define the signal at point Z
z_out = basebandports_to_pZ(Rx_baseband, target_found);
z = abs(z_out);


%% Find the largest correlator
% A method suggested in Class 9

% Define PN code
pn = [-1 -1 -1 1 1 -1 1];

% Sequence: a repetition of the PN code
pn_seq = [];
% Nc is equal to 7 and M is equal 199 (defined above)
for i=1:1:(Nc+1)*(M+1)
        pn_seq = [pn_seq pn];
end

% Find the largest correlator
c_in = 0;
largest_cor = 0;
k = 0;
for i=0:1:Nc
    if abs(xcorr(z,circshift(pn_seq,i),0)) > c_in
        c_in = abs(xcorr(z,circshift(pn_seq,i),0));
        largest_cor = c_in;
        k = i;
    end
end

% Define the sequence by using the largest correlator
sequence = circshift(pn_seq,k);

%% Define the Rx correlator: as a Matched Filter function
% Use of element-wise multiplication and summation

% One method is the convolution between signal z and PN code
% corrRx = conv(z,pn);
% However, this method might give wrong estimations

corrRx = [];

% Full correlation done manually as suggested in Class 9
for i=0:1:((Nc+1)*(M+1)-1)
    
    % Element-wise multiplication
    %window = sequence(i+1:i+7).*z(i+1:i+7);
    window = sequence(i+1:i+Nc).*z(i*Nc+1:i*Nc+Nc);
    
    % Summation
    sum_corr = sum(window);
    corrRx = [corrRx sum_corr];
end

% Define bin by bin for one PRI
bin = length(corrRx)/(Nc+1);
onePRI = zeros(1, bin);

for i=0:1:Nc
    RxPRI = corrRx(i*bin+1:i*bin+bin);
    onePRI = onePRI + RxPRI;
end

% Scaling: meters
N = 0:1:length(z)-1;
range = N*c*Tc/2;

% Averaging
rangeRx = range(1:Nc:end);
avgRx = onePRI./(Nc+1);
avgRange = rangeRx(1:bin);

%% Estimate the parameters of the target

Rx = abs(avgRx);
Amp = 0;

if targets == 1
    for i=1:1:M+1
        % Use of threshold to determine the peak
        if max(Rx(i))> threshold
            threshold = max(Rx(i));
            Amp = threshold;
            range_index = i;
        end
    end
end


if targets == 2
    peaks = 0;
    for i=1:1:M+1
        % Use of threshold to determine the peak
        if max(Rx(i))> threshold
            threshold  = max(Rx(i));
            Amp = threshold;
            range_index = i;
            peaks = [peaks Amp];
        end
    end
    Amp = sum(peaks)/(length(peaks));
end


% Estimate range at which the target occurs
R = avgRange(range_index);
techo = 2*R/c;

% Estimate the RCS value of the target
% Use of equations by using the method
% suggested in Class 8: Power mode
% Magnitude of the peak for one PRI

RCS = ((Amp/(45*45*1000))^2)*(4*pi)^3*(R^4)/(Ptx*(lambda^2));

%% 

% UNCOMMENT TO PLOT FIGURES

% figure();
% semilogy(z);
% title(['Signal at point Z due to Target-', num2str(targets)]);
% xlabel('Time Index (Tc)');
% ylabel('Magnitude (V)');
% xlim([0 11200]);
% ylim([1e-7 1e-3]);
% hold on;
% yline(v_thres,'--','Threshold', 'LineWidth', 1.5 ,'LabelHorizontalAlignment', 'Left');
% hold on;
% for i=1:1:(Nc+1)
%     xline(1400*i, 'red', ['PRI = ', num2str(i)], 'LineWidth', 1,'LabelHorizontalAlignment', 'Left');
% end
% set(gca, 'Fontsize', 14);
% 
% figure();
% semilogy(avgRange, abs(avgRx));
% title(['Range Detection for Target-', num2str(targets)]);
% xlabel('Range (m)');
% ylabel('Magnitude (V)');
% hold on;
% yline(v_thres,'--','Threshold', 'LineWidth', 1);
% set(gca, 'Fontsize', 14);

fprintf('Target: %d \n', targets);
fprintf('Target direction: %d degrees \n',target_found);
fprintf('Range: %d m',R);
fprintf(' = %d m \n',round(R));
fprintf('RCS average: %d m^2 ', RCS);
fprintf(' = %d m^2 \n',round(RCS));
fprintf('Range index: %dth\n',range_index);
fprintf('\n');

end