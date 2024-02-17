function show_cartesian(input,rr,r_lim,cr_lim,Nr,beam,caxis_range)
% * Show generated polar-format SLC image in Cartesian system
% * Inputs:
%  -input: magnitude image
%  -rr: range resolution in m
%  -r_lim: range axis limit in m. format: [r_min, r_max]
%  -cr_lim: cross-range axis limit in m. format: [cr_min, cr_max]
%  -Nr: number of samples in range direction.
%  -beam: antenna's cross-range beamwidth coverage in Radian.
%  -caxis_range: colormap limits

% Implemented by Benyamin Hosseiny
% Update: 2023-11-22

set(gcf,'renderer','zbuffer')

[fft_r,fft_cr] = size(input);
minRangeBinKeep = max(1,1+floor( (r_lim(1)/rr)*(fft_r/Nr) ));
maxRangeBinKeep = min(fft_r,1+floor( (r_lim(2)/rr)*(fft_r/Nr) ));

theta = -linspace( -rad2deg(beam)/2,rad2deg(beam)/2,fft_cr );
r_tick = linspace( rr,rr*Nr,fft_r )-rr;

[R_mat, theta_mat] = meshgrid(r_tick(minRangeBinKeep:maxRangeBinKeep), theta);
x_axis = R_mat.*cosd(theta_mat);
y_axis = R_mat.*sind(theta_mat);

input_clipped = flipud( input(minRangeBinKeep:maxRangeBinKeep,:)' );

if nargin<7
    caxis_range = [nanmin(input_clipped(:)), nanmax(input_clipped(:))];
end
if length(caxis_range)==1
    caxis_range = nanmax(input_clipped(:))+[-abs(caxis_range), 0];
end

fig = surf( y_axis, x_axis, input_clipped,'EdgeColor','none' ); view(2);
% fig = surf( y_axis, x_axis, input_clipped,'EdgeColor','k' ); view(2);
set( fig, 'AlphaData',~isnan(input_clipped) ) 
if caxis_range(1) ~= caxis_range(2)
	caxis( caxis_range );
end
axis('tight');axis('equal'); xlabel('X (m)'); ylabel('Y (m)'); 
xlim(cr_lim);ylim(r_lim); grid('on')
colormap('jet');%colormap('gray');
alpha 1

end