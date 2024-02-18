% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. "Structural displacement 
% monitoring using ground-based synthetic aperture radar." International Journal of 
% Applied Earth Observation and Geoinformation (JAG) 116 (2023): 103144.
% https://doi.org/10.1016/j.jag.2022.103144


clear;clc;close all
addpath('./src')
c = physconst('LightSpeed');
%% SYSTEM PARAMETERS
f       = 77e9;
lambda  = c/f;
fs      = 4.267e6;
Ts      = 60e-6;
bw      = 1e9;
snr     = 5; % (dB) Signal's SNR for adding white noise
Nr      = round(Ts*fs); % Range samples
prt     = 0.02;%2*Ts; % Data acquisition rate
prf     = 1/prt;
time    = 40;%.2; % Data acquisition time (sec)
V_radar0 = .02;            % Platform velocity (m/sec)
L_sar    = V_radar0*time;    % SAR Length (m)
Na       = 1+round(time*prf); % Azimuth samples

samples_time = round(time*prf);
r_freq  = 1/time; % frequency resolution!

%% Sensor positions:
X_rad = linspace(-L_sar/2, L_sar/2, Na);
Y_rad = 0*ones(size(X_rad));
Z_rad = 0*ones(size(X_rad));
% -------------------------------------------------------------------------
% % Stepper fluctuation function:
t_ax = linspace(0,time,Na);
dY_rad(1,:) = 0*(1e-3)*( exp(-t_ax).*cos(2*pi*(t_ax*2))+0.1*(rand(1,length(t_ax))-0.5) ); % Damped sine wave
% for ep_i = 1:epochs
%     dX_rad(ep_i,:) = ((-1)^ep_i)*(3e-3)*( exp(-t_ax).*cos(2*pi*(t_ax*2))+0.1*(rand(1,length(t_ax))-0.5) ); % Damped sine wave
% end
dX_rad = 0*repmat(t_ax,1,1);
dZ_rad = 0*repmat(t_ax,1,1);
% figure('Position', [20 100 1500 500]); sgtitle("Rail's fluctuations")
% subplot(1,3,1);plot(dX_rad(1,:));title('X');ylabel('Phase (rad)')
% subplot(1,3,2);plot(dY_rad(1,:));title('Y');ylabel('Phase (rad)')
% subplot(1,3,3);plot(dZ_rad(1,:));title('Z');ylabel('Phase (rad)')

% % Sensor positions, fluctuations included:
X_rad = X_rad+dX_rad; Y_rad = Y_rad+dY_rad; Z_rad = Z_rad+dZ_rad;

%% Target position:
X_tar = [0.02]; Y_tar = [2.1]; Z_tar = [-0.4];
% % Target Range with sensor's locations on the synthetic aperture:
R_tar_ref = sqrt( (X_rad-X_tar).^2 + (Y_rad-Y_tar).^2 + (Z_rad-Z_tar).^2  );
% % Adding 10e-6 gaussian noise:
R_tar = R_tar_ref+( 0*(10e-6)*2*(rand(size(R_tar_ref))-0.5) ); 
% % Target's backscattering amplitude:
amp_sig = 1; 

%% Target's displacement behaviour:
% Target's maximum displacement:
amp_defo_tar = 0.01; % (m) 
% Vibration frequency:
freq_tar     = .2;
% Sine displacement:
dR           = displacement_model_sin_TS(amp_defo_tar,freq_tar,prt,time)';
% Triangular displacement:
% period=2;
% dR           = displacement_model_triangle_TS(amp_defo_tar,period, prt,time)';

%% Target's Signal:
signal_Tar = signal_model_TS(amp_sig, Ts, fs, lambda, bw, R_tar, dR, snr);

%% Clutter:
N_clutter       = 20; % Number of clutters
sum_amp_clutter = 1.25*amp_sig;
amp_clutter     = (sum_amp_clutter/N_clutter)*ones(N_clutter,1); 
% amp_clutter     = 0.2*rand(N_clutter,1); 
R_clutter    = R_tar; % Same range with the target

% % Clutter's maximum displacement (non-zero if clutter is moving) :
amp_defo_clutter = 0*amp_defo_tar/30; % (m)
% % Clutter's Vibration frequency:
freq_clutter     = 5*freq_tar; 
% % Clutter's displacement:
dR_clutter       = displacement_model_sin_TS(amp_defo_clutter,freq_clutter,prt,time);

%% Clutter signal:
signal_clutter = signal_model_TS(amp_clutter, Ts, fs, lambda, bw, R_clutter, dR_clutter, snr);

% -------------------------------------------------------------------------
%% Raw data 2D array r*az [total contributions of clutters and target]:
cube_mimo_mono = signal_Tar'+signal_clutter';
cube_mimo_mono = awgn(cube_mimo_mono,snr);
%% ========================================================================
%% ========================================================================
%% ========================================================================
%% ========================================================================
