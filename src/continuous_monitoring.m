function [dR_hat_clutter_disturb_removed,dR_hat_clutter_removed,dR_hat_Rskew_removed,dR_hat_observed] = continuous_monitoring(rc_tar,lambda,prf,dstrb_rmv_flag,export_flag)
% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. "Structural displacement 
% monitoring using ground-based synthetic aperture radar." International Journal of 
% Applied Earth Observation and Geoinformation (JAG) 116 (2023): 103144.
% https://doi.org/10.1016/j.jag.2022.103144

% rc_tar: Detected target's range compressed signal during sensor's movement
% lambda: Signal's wavelenght
% dstrb_rmv_flag: Disturbance filtering flag
rc_tar=rc_tar(:);
% % Target's phase history
phi_rc_tar = angle(rc_tar); 
dR_hat_observed=(lambda/4/pi)*cumsum([0;wrapToPi( diff(phi_rc_tar) )]);
% -------------------------------------------------------------------------
%% 1-REMOVE Range skew (with reference phase or detrending):
% Estimated locations:
% % % X_tar_hat = r_x*round(X_tar/r_x); Y_tar_hat = r_y*round(Y_tar/r_y); Z_tar_hat = r_z*round(Z_tar/r_z);
% % % R_tar_hat = sqrt( (X_radar-X_tar_hat).^2 + (Y_radar-Y_tar_hat).^2 + (Z_radar-Z_tar_hat).^2  );% % Target Range with sensor's locations on the synthetic aperture
% % % phi_Rskew_tar = wrapToPi(4*pi*R_tar_hat/lambda)'; % Reference phase doesn't always work because that target's 3D location is limited to the digitized resolution cells
% % % phi_rc_tar_Rskew_removed = wrapToPi( phi_rc_tar-phi_Rskew_tar ); % Range skew compensated

% OR
% Detrending:
% phi_rc_tar_Rskew_removed = detrend( unwrap(phi_rc_tar),2 );
p=polyfit(linspace(1,length(phi_rc_tar),length(phi_rc_tar))',cumsum([0;wrapToPi( diff(phi_rc_tar) )]),2);
phi_rc_tar_Rskew_removed=cumsum([0;wrapToPi( diff(phi_rc_tar) )])-polyval(p,linspace(1,length(phi_rc_tar),length(phi_rc_tar))');
%
rc_tar_Rskew_removed = abs(rc_tar).*(cos(phi_rc_tar_Rskew_removed)+1i*sin(phi_rc_tar_Rskew_removed));
dR_hat_Rskew_removed = [0;wrapToPi(diff(angle(rc_tar_Rskew_removed)))];
dR_hat_Rskew_removed = (lambda/4/pi)*dR_hat_Rskew_removed;
dR_hat_Rskew_removed = cumsum(dR_hat_Rskew_removed);
% -------------------------------------------------------------------------
%% 2-REMOVE CLUTTERS:
% Circle fitting:
% circlefit_Par = CircleFitByPratt( [real(rc_tar),imag(rc_tar)] );
circlefit_Par = CircleFitByPratt( [real(rc_tar_Rskew_removed),imag(rc_tar_Rskew_removed)] );

cr_ampthresh = 0; % No threshold (suitable for long range or simulations)
cr_ampthresh = max(max(abs(rc_tar_Rskew_removed)))/1.5; % Suitable for close range of real experiments


idxabs=abs(rc_tar_Rskew_removed) >= cr_ampthresh;
rc_tar_Rskew_removed2=rc_tar_Rskew_removed(idxabs);
circlefit_Par = CircleFitByPratt( [real(rc_tar_Rskew_removed2),imag(rc_tar_Rskew_removed2)] );


x_c = circlefit_Par(1);
y_c = circlefit_Par(2);
radius_c = circlefit_Par(3);
% [x_c,y_c,radius_c,a] = circfit( real(rc(:,peak_c)),imag(rc(:,peak_c)) );
rc_clutter_removed = rc_tar_Rskew_removed-(x_c+1i*y_c);

dR_hat_clutter_removed = wrapToPi([0;diff(angle(rc_clutter_removed))]);
dR_hat_clutter_removed = (lambda/4/pi)*dR_hat_clutter_removed;
dR_hat_clutter_removed = cumsum(dR_hat_clutter_removed);
% -------------------------------------------------------------------------
%% 3-Filter Disturbances:
%[it is not the disturbances introduced by A_2021_Advancing...]
if dstrb_rmv_flag
    ff = fft(dR_hat_clutter_removed);
    xfig = 50; yfig = 0;   % Screen position
    widthfig  = 1000; heightfig = 500; % Width & Height of figure (by default in pixels)
    figure('Position', [xfig yfig widthfig heightfig]);
    
    freq_ax = linspace(0,prf,size(dR_hat_clutter_removed,1))';
    ff_dB = 10*log10( abs(ff(1:200)).^2 );
    % subplot(2,1,1); plot( dR_hat_clutter_removed, 'Color',[.5,1,.5], 'LineWidth', 1.5 )
    % subplot(2,1,1);
    plot( freq_ax(1:200), ff_dB, 'b', 'LineWidth', 2.5 ); hold on
    %max(ff_dB)-mean(ff_dB);
    
%     freq_idx = find( ff_dB==max(ff_dB) );
    [freq_sorted,freq_idx]=sort(ff_dB-max(ff_dB),'descend');
%     if any(freq_idx(1:2) == 1) % if first index exists in top 2
    if any(freq_idx(1) == 1) % if first index exists in top 2
        freq_idx=freq_idx(2);
    else
        freq_idx=freq_idx(1);
    end
    % Frequency domain filtering:
    ff(freq_idx+2:end-freq_idx-2)=0;
    if freq_idx~=1
        ff(1:max(freq_idx-3,1))=0;
        ff(end+3-max(freq_idx,3):end)=0;
    end
    % ff(freq_tar/r_freq+1+2:end-freq_tar/r_freq-1-2)=0;
    % ff(1:freq_tar/r_freq+1-3)=0;
    % ff(end-freq_tar/r_freq-1+3:end)=0;
    
    ff2=ifft(ff);
    dR_hat_clutter_disturb_removed = real(ff2);
    dR_hat_clutter_disturb_removed = dR_hat_clutter_disturb_removed-dR_hat_clutter_disturb_removed(1); % To make sure it starts from zero!
    ff_dB = 10*log10( abs(ff(1:200)).^2 );
    ff_dB(ff_dB==-inf)=-100;
    % subplot(2,1,2);
    plot( freq_ax(1:200), ff_dB, 'r', 'LineWidth', 1.1 );
    title('Displacment frequencies');xlabel('Frequency (Hz)');ylabel('Power (dB)')
    legend( 'Original', 'Filtered' );grid('on')
    if export_flag
        print(gcf,[export_directory 'Frequency response.jpg'],'-djpeg','-r600');
    end
else
    dR_hat_clutter_disturb_removed = [];
end
%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Phasor representation
% % Phasor Representation
figure('Position', [50 10 1500 500]); sgtitle('Phasor Representation')
% 1
subplot(1,4,1);scatter( real(rc_tar(idxabs)),imag(rc_tar(idxabs)),10, 'filled', 'MarkerFaceColor',[.5 .5 .5] );hold on; title('Observed (raw)')
% draw_circle(x_c,y_c,radius_c);hold on;
% scatter(x_c,y_c,'x','r','LineWidth',1.5);title('Observed signal')
error1 = mean( abs( abs(rc_tar)-radius_c ) )/(max(abs( rc_tar ))-min(abs(rc_tar)));

% 2
subplot(1,4,2);scatter( real(rc_tar_Rskew_removed(idxabs)),imag(rc_tar_Rskew_removed(idxabs)),10, 'filled', 'MarkerFaceColor',[.5 .5 .5] );hold on
draw_circle(x_c,y_c,radius_c);
scatter(x_c,y_c,'x','r','LineWidth',1.5); title('After correcting the range skew'); 
error2 = mean( abs(abs(rc_tar_Rskew_removed)-radius_c ) )/(max(abs( rc_tar_Rskew_removed ))-min(abs(rc_tar_Rskew_removed)));

% 3
subplot(1,4,3)
% rc_tar_clutter_removed = (rc_tar_Rskew_removed-x_c)+1i*(rc_tar_Rskew_removed-y_c);
rc_tar_clutter_removed = rc_tar_Rskew_removed-(x_c+1i*y_c);
scatter( real(rc_tar_clutter_removed(idxabs)), imag(rc_tar_clutter_removed(idxabs)),10, 'filled', 'MarkerFaceColor',[.5 .5 .5] );hold on
draw_circle(0,0,radius_c);
scatter(0,0,'x','r','LineWidth',1.5); title('After clutter reduction'); 
error3 = mean( abs(abs(rc_tar_clutter_removed)-radius_c) )/(max(abs(rc_tar_clutter_removed))-min(abs(rc_tar_clutter_removed)));
% error3 = mean( abs(abs(rc_tar_clutter_removed)-radius_c) )/mean(abs(rc_tar_clutter_removed));
% error3 = mean( abs(rc_tar_clutter_removed).^2)-radius_c^2;

% 4
% Reverse operation
zz0=diff(dR_hat_clutter_disturb_removed);zz1=zz0*4*pi/lambda;zz2=cumsum([0;zz1]);zz3=wrapToPi(zz2); % Cleaned phase of signal (Only Displacement component exists)
rc_tar_clutter_disturb_removed = abs(rc_tar_clutter_removed(1:end)).*(cos(zz3)+1i*sin(zz3));
circlefit4 = CircleFitByPratt( [real(rc_tar_clutter_disturb_removed),imag(rc_tar_clutter_disturb_removed)] );
x_c4= circlefit4(1); y_c4 = circlefit4(2); radius_c4 = circlefit4(3);
error4 = mean( abs(abs(rc_tar_clutter_disturb_removed)-radius_c4) )/(max(abs(rc_tar_clutter_disturb_removed))-min(abs(rc_tar_clutter_disturb_removed)));
% error4 = mean( abs(abs(rc_tar_clutter_disturb_removed)-radius_c4) )/mean(abs(rc_tar_clutter_disturb_removed));

subplot(1,4,4)
scatter( real(rc_tar_clutter_disturb_removed(idxabs)),imag(rc_tar_clutter_disturb_removed(idxabs)),10, 'filled', 'MarkerFaceColor',[.5 .5 .5] );hold on
title('After filtering'); 
% draw_circle(x_c4,y_c4,radius_c4)
% scatter(x_c4,y_c4,'x','r','LineWidth',1.5);title('After clutter and disturbance removal')

for ii = 1:4
    subplot(1,4,ii)
    box on; yline(0,'k--'); xline(0,'k--');
    axis('equal');grid on; xlim([-1.6,1.6]);ylim([-1.6,1.6]);xlabel('Real');ylabel('Imaginary')
end
if export_flag
    print(gcf,[export_directory 'Phasor representation.jpg'],'-djpeg','-r600');
end


function draw_circle(x,y,r)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit, 'r--', 'LineWidth',2);
end

end