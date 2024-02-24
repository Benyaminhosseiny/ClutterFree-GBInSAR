function signal_TS = signal_model_TS(Amp, Tp, fs, lambda, bw, R_tar, dR, snr)
% % % Used for simulating Continuous monitoring
% [1] A_2017_Static clutter elimination for frequency‐modulated continuous‐wave radar displacement
% [2] A_2019_An Effective Accuracy Evaluation Method for LFMCW Radar Displacement Monitoring with Phasor Statistical Analysis

% Amp   : Signal amplitude
% Tp    : Pulse(Sweep) time (sec)
% fs    : Sampling frequency
% lambda: central wavelength
% bw    : Bandwidth
% R_tar : Target's nominal range
% dR    : target's displacement over the time
% snr   : white noise snr (dB)

% -------------------------------------------------------------------------
c       = physconst('LightSpeed');
Nr      = round(Tp*fs);
t       = linspace(0,Tp,Nr);
fc      = c/lambda;
f0      = fc-bw/2;
lambda0 = c/f0;
fb_tar  = 2*bw*R_tar/c/Tp;

N_tar = length(Amp); % Number of Targets
N_ts  = length(dR); % Number of Time-series samples (Here also = Na)
for ii  = 1:N_ts
    signal_TS0 = 0;
    for jj = 1:N_tar
        signal_TS0 = signal_TS0+Amp(jj)*exp( -1i*(2*pi*fb_tar(ii)*t + 4*pi*R_tar(ii)/lambda0 + 4*pi*dR(ii)/lambda) ); % [1] Eq 1
        signal_TS(ii,:) = signal_TS0;
%         signal_TS(ii,:) = awgn( signal_TS0,snr ); % add Gaussian noise
    end
end
end