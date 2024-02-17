function ps_mask = PS(SLC_TS,thresh_adi,thresh_mag,coh_TS)
% PS detection:
% SLC_TS: r*c*t
% thresh_mag from max: dB

if nargin<3
    thresh_mag = 15; %dB
end
if nargin<2
    thresh_adi = 0.1;
end

% ADI
adi = ADI(SLC_TS);
% CDI
if nargin==4
    cdi = CDI(coh_TS);
else
    cdi=1;
end
% DI
DI = adi.*cdi;
% Mag
slc_abs = mean( 20*log10(abs(SLC_TS)),3 ); % dB
%
thresh_mag = max(max(slc_abs))-thresh_mag; %dB

mask_DI = DI; mask_DI(DI>thresh_adi) = 0; mask_DI(DI<=thresh_adi) = 1;

mask_mag = slc_abs; mask_mag(slc_abs<thresh_mag) = 0;mask_mag(slc_abs>=thresh_mag) = 1;

ps_mask = mask_DI.*mask_mag;

end