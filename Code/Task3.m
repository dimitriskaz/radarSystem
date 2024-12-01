%% Task 3: 1st scan - no targets

clc;
clear;
close all

%% Part a: Generate the noise snapshots at the baseband ports

%Zero targets are used

%Set a random theta_angle
theta_steer_angle = 30;

%Define signal at the baseband ports of Tx
Tx_baseband = pA_to_basebandports(theta_steer_angle);

%Define signal at the baseband ports of Rx
Rx_baseband = backscatterdata(Tx_baseband, 0); 

%% Part b: Plot the magnitude (Volts) of noise snapshots for one Dwell-time

z_out = basebandports_to_pZ(Rx_baseband,theta_steer_angle);
z = abs(z_out);

figure()
plot(z);
xlim([0 11200]);
ylim([0 2e-5]);
title('Noise for one Dwell-time');
xlabel('Snapshots');
ylabel('Magnitude (Volts)');
set(gca, 'Fontsize', 14);

%% Part c: Estimate and plot the pdf of the noise data samples

figure()
histogram(z,'Normalization','pdf');
title('PDF of the noise data samples');
xlabel('Magnitude (Volts)');
ylabel('PDF');
set(gca, 'Fontsize', 14);

%% Part d: Estimate the noise power at point Z

noise_sq = z.^2;
noise_power = mean(noise_sq);
