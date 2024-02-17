function adi = ADI(SLC_TS)
% SLC_TS: r*c*t
slc_ts_std = std( abs(SLC_TS),0,3 );
slc_ts_mean = mean( abs(SLC_TS),3 );
adi = slc_ts_std./slc_ts_mean;
end