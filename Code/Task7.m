%% Task 7: Radar data multi-target detection and parameter estimation

clc;
clear;
close all

%% Define parameters required

Tc = 28e-9;
Fc = 15e9;
c = 3e8;
lambda = c/Fc;
Nc = 7;
M = 199;
Ptx = 1;

%% Seach for the direction of the target (if any)

% NOTE: Download the BackscatterData.mat from the dk3617 file
% because it is not included in my .zip file
% the BackscatterData.mat file was a very large file and my 
% dk3617.zip file could not be submitted by including it
load BackscatterData.mat
    
data = BackscatterData;

% Convert the data cell to matrices
z_out = [];
for i = 1:1:121
    Rx_baseband = cell2mat(data([i]));
    theta_steer_angle = i+29;
    z_out(i,:) = basebandports_to_pZ(Rx_baseband,theta_steer_angle);
end

% Z is a 121x11,200 matrix, which means that each row (1x11,200)
% has stored the data for each angle
% e.g. z(1,:) = 1x11,200 matrix data for angle = 30
%      z(2,:) = 1x11,200 matrix data for angle = 31 etc.

z = abs(z_out); 

% Set the threshold
max_val = 1e-5;
threshold = max_val;

index = 0;

for i = 1:1:121
    % Check if a target is detected
    maximum = max(z(i,:));
    if maximum > max_val
        max_val = maximum;
        index = i;
    end
end

% The azimuth angle takes values from 30 to 150:
target_found = index + 29;
% Define the signal according to the direction found
Rx_baseband = cell2mat(data([index]));
z_out = basebandports_to_pZ(Rx_baseband,target_found);
z = abs(z_out);
%plot(z);


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
    window = sequence(i+1:i+Nc).*z(i*Nc+1:i*Nc+Nc);
    %window = sequence(i+1:i+7).*z(i+1:i+7);
    
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

%% Estimate parameters of the target

Rx = abs(avgRx);
figure();
plot(Rx);
title('Rx correlator');
xlabel('Number of bins');
ylabel('Magnitude (V)');
set(gca, 'Fontsize', 14);

num= 1;
Amp = [];
range_index = [];

for i=1:1:M+1
    % Use of threshold to determine the peaks
    if max(Rx(i))> 1e-5
        Amp(num) = max(Rx(i));
        range_index(num) = i;
        num = num +1;   
    end
end

% Target 1
Amp1 = Amp(2);
range_index1 = range_index(2);
R1 = avgRange(range_index1);
techo1 = 2*R1/c;
RCS1 = ((Amp1/(45*45*1000))^2)*(4*pi)^3*(R1^4)/(Ptx*(lambda^2));

% Target 2
Amp2 = Amp(3);
range_index2 = range_index(3);
R2 = avgRange(range_index2);
techo2 = 2*R2/c;
RCS2 = ((Amp2/(45*45*1000))^2)*(4*pi)^3*(R2^4)/(Ptx*(lambda^2));

% Target 3
Amp3 = Amp(4);
range_index3 = range_index(4);
R3 = avgRange(range_index3);
techo3 = 2*R3/c;
RCS3 = ((Amp3/(45*45*1000))^2)*(4*pi)^3*(R3^4)/(Ptx*(lambda^2));

% Target 4
Amp4 = Amp(6);
range_index4 = range_index(6);
R4 = avgRange(range_index4);
techo4 = 2*R4/c;
RCS4 = ((Amp4/(45*45*1000))^2)*(4*pi)^3*(R4^4)/(Ptx*(lambda^2));

% %UNCOMMENT TO PLOT FIGURES

% figure();
% semilogy(z);
% title('Signal at point Z');
% xlabel('Time Index (Tc)');
% ylabel('Magnitude (V)');
% xlim([0 11200]);
% ylim([1e-7 1e-3]);
% hold on;
% yline(threshold,'--','Threshold', 'LineWidth', 1.5 ,'LabelHorizontalAlignment', 'Left');
% hold on;
% for i=1:1:(Nc+1)
%     xline(1400*i, 'red', ['PRI = ', num2str(i)], 'LineWidth', 1.5,'LabelHorizontalAlignment', 'Left');
% end
% set(gca, 'Fontsize', 14);
% 
% figure();
% semilogy(avgRange, abs(avgRx));
% title('Range Detection');
% xlabel('Range (m)');
% ylabel('Magnitude (V)');
% hold on;
% yline(threshold,'--','Threshold', 'LineWidth', 1.5 ,'LabelHorizontalAlignment', 'Left');
% set(gca, 'Fontsize', 14);

fprintf('-- Target 1 -- \n');
fprintf('Target direction: %d degrees \n',target_found);
fprintf('Range: %d m',R1);
fprintf(' = %d m \n',round(R1));
fprintf('RCS average: %d m^2 \n', RCS1);
fprintf('Range index: %dth\n',range_index1);
fprintf('\n');


fprintf('-- Target 2 -- \n');
fprintf('Target direction: %d degrees \n',target_found);
fprintf('Range: %d m',R2);
fprintf(' = %d m \n',round(R2));
fprintf('RCS average: %d m^2 ', RCS2);
fprintf(' = %d m^2 \n',round(RCS2));
fprintf('Range index: %dth\n',range_index2);
fprintf('\n');


fprintf('-- Target 3 -- \n');
fprintf('Target direction: %d degrees \n',target_found);
fprintf('Range: %d m',R3);
fprintf(' = %d m \n',round(R3));
fprintf('RCS average: %d m^2 ', RCS3);
fprintf(' = %d m^2 \n',round(RCS3));
fprintf('Range index: %dth\n',range_index3);
fprintf('\n');

fprintf('-- Target 4 -- \n');
fprintf('Target direction: %d degrees \n',target_found);
fprintf('Range: %d m',R4);
fprintf(' = %d m \n',round(R4));
fprintf('RCS average: %d m^2', RCS4);
fprintf(' = %d m^2 \n',round(RCS4));
fprintf('Range index: %dth\n',range_index4);
fprintf('\n');
