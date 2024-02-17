% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. "Structural displacement 
% monitoring using ground-based synthetic aperture radar." International Journal of 
% Applied Earth Observation and Geoinformation (JAG) 116 (2023): 103144.
% https://doi.org/10.1016/j.jag.2022.103144

% * This code is developed for reducing the effect of static clutter from the TS of radar cube data [r x az x time]
% * run "ClutterR_Part2_processing.m" after running this code!

% * First, run "ClutterR_Part2_processing.m"
% * Or
% * import your own "cube_3dTS" (a 3D array of TS of raw cube data [r x az x time])

% clear;clc;close all
c = physconst('LightSpeed');
addpath('./src')

%% Input setting:
% % SAR Processing:
r_fft  = size(cube_3dTS,1)/2;
az_fft = 64;
el_fft = 1;

% % InSAR Processing:
adi_thresh = .5;
amp_thresh = 1; % dB from the max

% =========================================================================

%% Known Parameters:
Nr = size(cube_3dTS,1);
Na = size(cube_3dTS,2);
Ne = size(cube_3dTS,3);
Nts= size(cube_3dTS,4); % Number of TS data
d_az = d_az;
d_el = d_el;

beam_az  = rad2deg(lambda/2/d_az);
beam_el  = rad2deg(lambda/2/d_el);
R_ax     = linspace(rr/2,rr*Nr-rr/2,r_fft)'; 
Theta_ax = -linspace(-beam_az/2,beam_az/2,az_fft);
Phi_ax   = -linspace(-beam_el/2,beam_el/2,el_fft);

% Radar antennas pos:
X_rad = d_az*[-floor(Na/2):floor(Na/2)];
Y_rad0 = 0;
Z_rad0 = d_el*[-floor(Ne/2):floor(Ne/2)];
nadir=-0;
rot_mat=[1 0 0; 0 cosd(nadir) -sind(nadir); 0 sind(nadir) cosd(nadir)];
z_shift=0;

Y_rad = Y_rad0*cosd(nadir)-Z_rad0*sind(nadir);
Z_rad = z_shift+Y_rad0*sind(nadir)+Z_rad0*cosd(nadir);

%% ========================================================================
%% ========================================================================
%% ========================================================================
%% ========================================================================

%% SAR Processing:
rcTS  = fft( cube_3dTS,r_fft,1 );
slc3dTS = fftshift( fft( rcTS,az_fft,2 ),2 );
slc3dTS = slc3dTS./max(abs(slc3dTS),[],'all');
slc3dTS = squeeze(slc3dTS);

slc_db = 10*log10(abs(slc3dTS(:,:,10))); 
slc_db = slc_db-max(slc_db,[],'all'); % Set maximum backscatter to zero
clear cube_3dTS

%% Plot 1 SLC:
% Define plotting characteristics:
az_lim = (beam_az/2)*[-1, 1];        % azimuth axis limit in polar format plot (Deg)
r_lim  = [0, 50];                    % Range axis limit in polar format plot (m)
cr_lim = [-20, 20];                  % Cross-Range axis (X) limit in Cartesian format plot (m)
caxis_range = [-20, 0]; % Colormap limit
% Show in Polar format:
figure("Position",[50,50,1500,700]);
subplot(1,2,1); show_polar(slc_db, rr, r_lim, az_lim, Nr, deg2rad(beam_az), caxis_range); colorbar()
title('SAR image in Polar Coordinate')

% Show in Cartesian format:
subplot(1,2,2); show_cartesian(slc_db, rr, r_lim, cr_lim, Nr, deg2rad(beam_az), caxis_range); colorbar()
title('SAR image in Cartesian Coordinate')
%% InSAR Processing:
% % %  ADI analysis:
adi = ADI(slc3dTS);
% % % PS detection:
ps_mask = PS(slc3dTS, adi_thresh, amp_thresh);%ps_mask(ps_mask==0)=NaN;

% % % point cloud:
[ps_idxr, ps_idxaz] = find(ps_mask==1);

%% TSInSAR:
%% Clutter reduction
for ii=1:length(ps_idxr)
    PS_pixelii        = slc3dTS( ps_idxr(ii),ps_idxaz(ii),: );
    PS_pixel(ii,:)    = PS_pixelii;
    [~,PS_ClutterR(ii,:)] = Clutter_reduction( PS_pixelii(:) );
    ps_Theta(ii,1)= Theta_ax(ps_idxaz(ii));
    ps_R(ii,1)    = R_ax(ps_idxr(ii));
    ps_adi(ii,1)  = adi(ps_idxr(ii),ps_idxaz(ii));
    ps_mag(ii,1)  = slc_db(ps_idxr(ii),ps_idxaz(ii));
end

PS_defo    = cumsum( wrapToPi( diff( angle(PS_pixel),1,2 ) ),2 )*lambda/4/pi;    % Observed deformation signal (before clutter reduction)
PS_defo_CR = cumsum( wrapToPi( diff( angle(PS_ClutterR),1,2 ) ),2 )*lambda/4/pi; % Corrected deformation signal (after clutter reduction)


%% Phasor representation
ii = 1;
[circle_par, PS_ClutterRii] = Clutter_reduction( PS_pixel(ii,:) );
% % Phasor Representation
figure('Position', [50 10 1500 1000]); sgtitle('Clutter Reduction')

subplot(3,3,3);scatter( real(PS_pixel(ii,:)), ...
                        imag(PS_pixel(ii,:)), 10, 'filled', 'MarkerFaceColor',[.5 .5 .5] );hold on; title("before clutter reduction")
draw_circle( circle_par(1), circle_par(2), circle_par(3) )
subplot(3,3,6);scatter( real(PS_ClutterR(ii,:)), ...
                        imag(PS_ClutterR(ii,:)), 10, 'filled', 'MarkerFaceColor',[.5 .5 .5] );hold on; title("after clutter reduction")
for pp = [3,6]
    subplot(3,3,pp)
    box on; yline(0,'k--'); xline(0,'k--');
    axis('equal');grid on; xlim([-1,1]);ylim([-1,1]);xlabel('Real');ylabel('Imaginary')
end                    
subplot(3,3,[1,2]); plot(PS_defo(ii,:)',    'b'); title("raw displacement signal (before clutter reduction"); 
ylim( max(abs([PS_defo_CR(:); PS_defo(:)]))*[-1,1] );ylabel('displacement [m]');xlabel('Time series samples')
subplot(3,3,[4,5]); plot(PS_defo_CR(ii,:)', 'b'); title("corrected displacement signal (after clutter reduction");
ylim( max(abs([PS_defo_CR(:); PS_defo(:)]))*[-1,1] );ylabel('displacement [m]');xlabel('Time series samples')
subplot(3,3,[7,8]); plot(dR_tar(:,1),       'b'); title("Reference");      
ylim( max(abs([PS_defo_CR(:); PS_defo(:)]))*[-1,1] );ylabel('displacement [m]');xlabel('Time series samples')     


% =========================================================================
function draw_circle(x,y,r)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit, 'r--', 'LineWidth',2);
end