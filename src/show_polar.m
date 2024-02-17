function show_polar(input,rr,r_lim,az_lim,Nr,beam,caxis_range)
% * Show generated polar-format SLC image in Cartesian system
% * Inputs:
%  -input: magnitude image
%  -rr: range resolution in m
%  -r_lim: range axis limit in m. format: [r_min, r_max]
%  -ar_lim: azimuth axis limit in degree. format: [az_min, az_max]
%  -Nr: number of samples in range direction.
%  -beam: antenna's cross-range beamwidth coverage in Radian.
%  -caxis_range: colormap limits

% Implemented by Benyamin Hosseiny
% Update: 2023-11-22

[fft_r, fft_cr] = size(input);
theta_tick = -linspace( -rad2deg(beam)/2,rad2deg(beam)/2,fft_cr );
r_tick = linspace( rr,rr*Nr,fft_r )-rr;

minRangeBinKeep = max(1,1+floor( (r_lim(1)/rr)*(fft_r/Nr) ));
maxRangeBinKeep = min(fft_r,0+floor( (r_lim(2)/rr)*(fft_r/Nr) ));
minThetaBinKeep = max(1,1+floor( (fft_cr/2)*(1+az_lim(1)/(rad2deg(beam)/2)) ));
maxThetaBinKeep = min(fft_cr,0+floor( (fft_cr/2)*(1+az_lim(2)/(rad2deg(beam)/2)) ));

input_clipped = fliplr( input(minRangeBinKeep:maxRangeBinKeep,minThetaBinKeep:maxThetaBinKeep) );
if nargin<7
    caxis_range = [min(min(input_clipped)), max(max(input_clipped))];
end
if length(caxis_range)==1
    caxis_range = max(input_clipped,[],'all')+[-abs(caxis_range), 0];
end

fig = imagesc(theta_tick(minThetaBinKeep:maxThetaBinKeep),r_tick(minRangeBinKeep:maxRangeBinKeep), input_clipped);set(gca,'YDir','normal'); 
set( fig, 'AlphaData',~isnan(input_clipped) )
colormap('jet');xlabel('\theta (deg)'); ylabel('R (m)'); grid('on'); 

if caxis_range(1) ~= caxis_range(2)
    caxis( caxis_range );%xlim(az_lim)
end

end