% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. "Structural displacement 
% monitoring using ground-based synthetic aperture radar." International Journal of 
% Applied Earth Observation and Geoinformation (JAG) 116 (2023): 103144.
% https://doi.org/10.1016/j.jag.2022.103144

% * This code simulates the process of reducing the effect of static clutter from the TS of radar cube data [r x az x time]
% * run "ClutterR_Part2_processing.m" after running this code!


clear;
clc;
close all
c = physconst('LightSpeed');
addpath('./src')
%% Parameters:
fc = 77e9; 
lambda = c/fc;
bw = .3e9;         % Signal bandwdith
rr = c/2/bw;       % Range resolution
T = 60e-6; 
cr = bw/T;
prt = 0.02;        % Data repetition time
time= 20;          % Observation time
snr = 0;
Nr = 256;          % Number of range samples
Na = 15;           % Number of azimuth samples
Ne = 1;            % Number of elevation samples
Nts= time/prt;           % Number of TS data

fs = Nr/T;
d_az = lambda/4;
d_el = 0;2*lambda/2;

%% Radar position:
X_rad0  = d_az*(0:Na-1)-d_az*Na/2; %d_az*[-floor(Na/2):floor(Na/2)];
Y_rad0  = zeros(Na, Ne);
Z_rad0  = d_el*(0:Ne-1)-d_el*Ne/2; %d_el*[-floor(Ne/2):floor(Ne/2)];

% Define rotation around X axis (nadir or incidence angle):
nadir   = -0;
rot_mat = [1 0 0; 0 cosd(nadir) -sind(nadir); 0 sind(nadir) cosd(nadir)];
z_shift = 0;

% Sensor's moving surface:
Y_rad = Y_rad0*cosd(nadir)-Z_rad0*sind(nadir);         % shape: Na*Ne
Z_rad = z_shift+Y_rad0*sind(nadir)+Z_rad0*cosd(nadir); % shape: Na*Ne
X_rad = X_rad0(:)*ones(1,Ne);                          % shape: Na*Ne

figure('Position', [500 500 800 500]); 
scatter(X_rad(:),Z_rad(:)); axis equal; title("Sensor's trajectory"); xlabel('Azimuth (m)'); ylabel('Elevation (m)'); grid on; box on
% scatter3(Y_rad(:),X_rad(:),Z_rad(:)); axis auto; title("Sensor's trajectory"); xlabel('Range (m)'); ylabel('Azimuth (m)'); zlabel('Elevatio (m)')
xlim(1.1*[min(X_rad),max(X_rad)])

% Reshape it to 3D cube:
X_rad= repmat( reshape(X_rad, 1,Na,Ne), Nr,1,1); % shape: Nr*Na*Ne
Y_rad= repmat( reshape(Y_rad, 1,Na,Ne), Nr,1,1); % shape: Nr*Na*Ne
Z_rad= repmat( reshape(Z_rad, 1,Na,Ne), Nr,1,1); % shape: Nr*Na*Ne

% Add stepper fluctuation function:
xrad_err=0;2*1e-5;
yrad_err=0;2*1e-5;
zrad_err=0;2*1e-5;

%% Target position:
% Polar
% R_tar=[2]; theta_tar=[20.2]; ph_tar=[-20.5];
% X_tar=R_tar.*sind(theta_tar).*cosd(ph_tar);
% Y_tar=R_tar.*cosd(theta_tar).*cosd(ph_tar);
% Z_tar=R_tar.*sind(ph_tar);

% % Cartesian:
X_tar=[ 0.1];
Y_tar=[10];
Z_tar=[ 0];
A_tar=[1]; % target backscattering amplitude

%% Clutter position:
N_clutter = 10;
% sum_amp_clutter = 15.25*A_tar;
% amp_clutter     = (sum_amp_clutter/N_clutter)*ones(1,N_clutter); 

X_tar= X_tar+[ 0,  0.1*ones(1, N_clutter).*(rand(1, N_clutter)-0.5)];
Y_tar= Y_tar+[ 0,  0.1*ones(1, N_clutter).*(rand(1, N_clutter)-0.5)];
Z_tar= Z_tar+[ 0,  0.1*ones(1, N_clutter).*(rand(1, N_clutter)-0.5)];
A_tar= A_tar+[ 0,  0.5+0.1*ones(1, N_clutter).*(rand(1, N_clutter)-0.5)];
% A_tar= [ A_tar,  amp_clutter];

R_tar    =sqrt(X_tar.^2+Y_tar.^2+Z_tar.^2);
theta_tar=asin(X_tar./R_tar);
ph_tar   =asin(Z_tar./R_tar);
N_tar=length(A_tar);

figure("Position", [100,100,1000,700]); 
subplot(1,2,1);
scatter(X_tar,Y_tar); xlim([-1.1*max(abs(X_tar)+.1), 1.1*max(abs(X_tar)+.1)]); ylim([-0.1, 1.1*max(abs(Y_tar)+.1)]); grid on; axis equal
title('Targets location'); xlabel('X (m)'); ylabel('Y (m)'); box on
subplot(1,2,2);
scatter(X_tar,Y_tar,20); xlim([-0.5+min(X_tar), 0.5+max(X_tar)]); ylim([-1+min(Y_tar), 1+max(Y_tar)]); grid on; axis equal
title('Targets location (Zoomed)'); xlabel('X (m)'); ylabel('Y (m)'); box on

%% Displacement per TS:
dX_tar = [1:Nts]'*zeros(1, N_tar);
dZ_tar = [1:Nts]'*zeros(1, N_tar);
dY_tar = [1:Nts]'*[-0.5e-3, zeros(1, N_tar-1)]; % Line of Sight

amp_defo_tar = 0.01; % (m) 
% Vibration frequency:
freq_tar = 0.2;
% Sine displacement:
dY_tar0 = displacement_model_sin_TS(amp_defo_tar,freq_tar,prt,time);
dY_tar  = [dY_tar0(2:end), [1:Nts]'*zeros(1, N_tar-1)];

dR_tar  = dY_tar;

t = linspace(0,T,Nr)'; % Pulse (sweep) time axis
t = repmat(t, 1, Na, Ne); % shape: Nr*Na*Ne

%% Raw cube SIGNAL:
for ts_ii = 1:Nts
    % Add stepper fluctuation:
    X_rad2 = X_rad + xrad_err*(rand(1, Na, Ne)-0.5);
    Y_rad2 = Y_rad + yrad_err*(rand(1, Na, Ne)-0.5);
    Z_rad2 = Z_rad + zrad_err*(rand(1, Na, Ne)-0.5);

    %Rail's starting point error (spe) in X axis (unsynchronization):
    spe_shift=1; % shift in the steps (e.g., 1 step)
    spe_X(:,ts_ii) = 0*(-1.^randi(2,Ne,1).*randi(spe_shift,Ne,1).*lambda/8); 
    X_rad2 = X_rad2+reshape(spe_X(:,ts_ii),1,1,Ne); % + Starting point error (mismatch)
    
    cube_3d_ii=0;
    for tar_ii = 1:N_tar
        X_tarii = X_tar(tar_ii)+dX_tar(ts_ii, tar_ii);
        Y_tarii = Y_tar(tar_ii)+dY_tar(ts_ii, tar_ii);
        Z_tarii = Z_tar(tar_ii)+dZ_tar(ts_ii, tar_ii);
        
% %         % Signal Model in Polar System:
%         R_tar_ii = R_tar(tar_ii)+dR_tar(ts_ii,tar_ii);
%         tau = 2*( R_tar_ii )/c;
%         cube_3d_ii = cube_3d_ii + exp( 1i*2*pi*( fc*tau + cr*t.*tau - cr*(tau.^2)/2 - ...
%                                                  2*X_rad2.*sind(theta_tar(tar_ii))/lambda - ...
%                                                  2*Z_rad2.*sind(ph_tar(tar_ii))/lambda ) );
        

        % Signal Model in Cartesian System:
        R_tar_ii = ...
            sqrt( ( X_tarii-X_rad2 ).^2 + ...
                  ( Y_tarii-Y_rad2 ).^2 + ...
                  ( Z_tarii-Z_rad2 ).^2 );

        tau = 2*( R_tar_ii )/c;   
        cube_3d_ii = cube_3d_ii + A_tar(tar_ii)*exp( 1i*2*pi*( fc*tau + cr*t.*tau - cr*(tau.^2)/2 ) );
        
    end
    cube_3dTS(:,:,:,ts_ii) = awgn( cube_3d_ii, snr );
end
