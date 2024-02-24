function [rc,slc] = create_slc(cube_mimo_mono,r_fft,az_fft)
% cube_mimo_mono: range_dir x cr_dir

rc  = fft( cube_mimo_mono,r_fft,1 );
slc = fftshift( fft( rc,az_fft,2 ),2 );
end