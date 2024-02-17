function [circle_par, sig_clutter_removed] = Clutter_reduction(sig_cmplx)
% Circle fitting:
sig_cmplx = sig_cmplx(:);
idxabs=abs(sig_cmplx)>max(max(abs(sig_cmplx)))/1.5;%>0;
sig_cmplx2=sig_cmplx(idxabs);
circlefit_Par = CircleFitByPratt( [real(sig_cmplx2),imag(sig_cmplx2)] );
x_c = circlefit_Par(1);
y_c = circlefit_Par(2);
radius_c = circlefit_Par(3);
% [x_c,y_c,radius_c,a] = circfit( real(rc(:,peak_c)),imag(rc(:,peak_c)) );
sig_clutter_removed = sig_cmplx-(x_c+1i*y_c);
circle_par = [x_c,y_c,radius_c];
end