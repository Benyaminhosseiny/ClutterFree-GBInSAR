% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. "Structural displacement 
% monitoring using ground-based synthetic aperture radar." International Journal of 
% Applied Earth Observation and Geoinformation (JAG) 116 (2023): 103144.
% https://doi.org/10.1016/j.jag.2022.103144
% * First, run "Continuous_dsplcmnt_Part1_rawdata.m"
% * Or
% * import your own "cube_mimo_mono" (a 2D array of single-track SAR data [r x az])

addpath('./src')
c = physconst('LightSpeed');

%% Known PARAMETERS
% f       = 77e9;
% Ts      = 60e-6;
% bw      = 1e9;
% prt     = 0.02;                  % Data acquisition rate (sec)
% time    = 40;                    % Data acquisition time (sec)
% V_radar0= .025;                  % Platform velocity (m/sec)
% L_sar    = V_radar0*time;        % SAR Length (m)
% Nr = size(cube_mimo_mono,1);     % Range samples
% Na = size(cube_mimo_mono,2);     % Azimuth samples
% lambda  = c/f;
% prf     = 1/prt;
% samples_time = round(time*prf);
% r_freq  = 1/time;                % frequency resolution!

%% Input setting:
r_fft=Nr*2;
az_fft=Na;
dstrb_rmv_flag = 1;
d_az     = prt*V_radar0;         % Antenna spacing (Azimuth sampling)


%% ========================================================================
%% ========================================================================
%% ========================================================================
%% ========================================================================

%% Target detection & Range Compression:
[rc,slc] = create_slc(cube_mimo_mono,r_fft,az_fft);
rc = rc./max(abs(rc),[],'all');

% Finding target's peak point
[peak_r,~] = find(rc==max(max(rc)));
% peak_r=50;26*2;54;51;%
rc_tar = rc(peak_r,:);

%% continuous monitoring:
export_flag    = 0;
[dR_hat_clutter_disturb_removed,dR_hat_clutter_removed,dR_hat_Rskew_removed,dR_hat_observed] = continuous_monitoring(rc_tar,lambda,prf,dstrb_rmv_flag,export_flag);

% Zero mean:
dR_hat_observed = dR_hat_observed-mean(dR_hat_observed);
dR_hat_Rskew_removed = dR_hat_Rskew_removed-mean(dR_hat_Rskew_removed);
dR_hat_clutter_removed = dR_hat_clutter_removed-mean(dR_hat_clutter_removed);
dR_hat_clutter_disturb_removed = dR_hat_clutter_disturb_removed-mean(dR_hat_clutter_disturb_removed);

%%
% figure
% subplot(2,2,1);imagesc(abs(rc)); title('range compressed signal')
% subplot(2,2,2);imagesc(abs(slc));title('slc')
% subplot(2,2,3);plot(abs(rc_tar)); title('target range compressed signal')

s1=50;1;150;
s2=1100;Na;950;
xax=prt*(0:s2-s1);
figure("Position", [50 50 1000 1000]);
subplot(5,1,1);plot(xax,dR_hat_observed(s1:s2)*1e3,       'b','LineWidth',1.5);         title('Displacement from observed signal');       xlabel('Time [s]');ylabel('[mm]')
subplot(5,1,2);plot(xax,dR_hat_Rskew_removed(s1:s2)*1e3,  'b','LineWidth',1.5);         title('Displacement after deskew');               xlabel('Time [s]');ylabel('[mm]');ylim([-15,15])
subplot(5,1,3);plot(xax,dR_hat_clutter_removed(s1:s2)*1e3,'b','LineWidth',1.5);         title('Displacement after clutter removal');      xlabel('Time [s]');ylabel('[mm]');ylim([-15,15])
if dstrb_rmv_flag
subplot(5,1,4);plot(xax,dR_hat_clutter_disturb_removed(s1:s2)*1e3,'b','LineWidth',1.5); title('Displacement after disturbance filtering');xlabel('Time [s]');ylabel('[mm]');ylim([-15,15])
end
if exist('dR','var')
    subplot(5,1,5);plot(xax,dR(s1:s2)*1e3,'b','LineWidth',1.5); title('Reference');xlabel('Time [s]');ylabel('[mm]');ylim([-15,15])
end

%% New figures:
zoom1 = s1+740;855;855;750;
zoom2 = s1+790;897;895;950;
figure('Position', [50 0 1500 1500]);
subplot(4,4,[1,2])
plot( t_ax(s1:s2), 1e3*dR_hat_observed(s1:s2), 'b', 'LineWidth', 2 );hold on
plot( t_ax(s1:s2), 1e3*dR(s1:s2), 'r' )
ylabel('[mm]');title('Radar captured displacement')
ylim(1.1*max(abs(1e3*[dR(s1:s2),dR_hat_observed(s1:s2)',dR_hat_Rskew_removed(s1:s2)',dR_hat_clutter_removed(s1:s2)',dR_hat_clutter_disturb_removed(s1:s2)']))*[-1,1])

subplot(4,4,[5,6])
plot( t_ax(s1:s2), 1e3*dR_hat_Rskew_removed(s1:s2), 'LineWidth', 2, 'Color',[.5,.5,.5] );hold on
plot( t_ax(s1:s2), 1e3*dR(s1:s2), 'r' )
ylabel('[mm]');title('Displacement after range deskew')
ylim(1.1*max(1e3*abs([dR(s1:s2),dR_hat_observed(s1:s2)',dR_hat_Rskew_removed(s1:s2)',dR_hat_clutter_removed(s1:s2)',dR_hat_clutter_disturb_removed(s1:s2)']))*[-1,1])

subplot(4,4,[9,10])
plot( t_ax(s1:s2), 1e3*dR_hat_clutter_removed(s1:s2), 'LineWidth', 2, 'Color',[.5,1,.5] );hold on
plot( t_ax(s1:s2), 1e3*dR(s1:s2), 'r' )
rectangle('Position',[ .99*t_ax(zoom1), 1e3*dR(zoom1), t_ax(zoom2)-.97*t_ax(zoom1), 1.7*1e3*(max(dR(s1:s2))-dR(zoom1)) ],'EdgeColor','b','LineStyle','--','LineWidth',1);
ylabel('[mm]');title('Displacement after clutter reduction')
ylim(1.1*max(abs(1e3*[dR(s1:s2),dR_hat_observed(s1:s2)',dR_hat_Rskew_removed(s1:s2)',dR_hat_clutter_removed(s1:s2)',dR_hat_clutter_disturb_removed(s1:s2)']))*[-1,1])

subplot(4,4,[13,14])
plot( t_ax(s1:s2), 1e3*dR_hat_clutter_disturb_removed(s1:s2), '--', 'LineWidth', 2, 'Color', [0, 0.0, 0.0] );hold on
plot( t_ax(s1:s2), 1e3*dR(s1:s2), 'r' );hold('on')
rectangle('Position',[ .99*t_ax(zoom1), 1e3*dR(zoom1), t_ax(zoom2)-.97*t_ax(zoom1), 1.7*1e3*(max(dR(s1:s2))-dR(zoom1)) ],'EdgeColor','b','LineStyle','--','LineWidth',1);
ylim(1.1*max(abs(1e3*[dR(s1:s2),dR_hat_observed(s1:s2)',dR_hat_Rskew_removed(s1:s2)',dR_hat_clutter_removed(s1:s2)',dR_hat_clutter_disturb_removed(s1:s2)']))*[-1,1])
ylabel('[mm]');xlabel('Time [s]');title('Displacement after disturbance filtering')

subplot(4,4,[3,4,7,8])
plot( t_ax(zoom1:zoom2), 1e3*dR(zoom1:zoom2), 'r', 'LineWidth', 3 );hold on
plot( t_ax(zoom1:zoom2), 1e3*dR_hat_observed(zoom1:zoom2), 'b', 'LineWidth', 2 );hold on
plot( t_ax(zoom1:zoom2), 1e3*dR_hat_Rskew_removed(zoom1:zoom2), 'LineWidth', 2, 'Color',[.5,.5,.5] );hold on
plot( t_ax(zoom1:zoom2), 1e3*dR_hat_clutter_removed(zoom1:zoom2), 'LineWidth', 2.5, 'Color',[.5,1,.5] );hold on
plot( t_ax(zoom1:zoom2), 1e3*dR_hat_clutter_disturb_removed(zoom1:zoom2), '--', 'LineWidth', 2, 'Color', [0, 0.0, 0.0] );grid('on')
ylim( 1e3*[0.9*min(dR(zoom1:zoom2)-1e-30),1.1*max(dR(zoom1:zoom2)+1e-30)] )
ylabel('Displacement [mm]');xlabel('Time [s]');title('Zoomed area')

h = legend('reference','observed raw displacement signal','after correcting range skew','after clutter removal', 'after filtering');
pos = get(h,'Position');
posx = 0.7;
posy = 0.2;
set(h,'Position',[posx posy pos(3) pos(4)]);

