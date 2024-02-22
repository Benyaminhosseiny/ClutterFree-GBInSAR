function dR_sin = displacement_model_sin_TS(amp,freq,prt,time,fi_shift)
% amp : displacement amplitude (m)
% freq: vibration frequency
% prt : data acquisition (sampling) interval
% time: data acquisition (sampling) total time
% fi_shift: starting phase shift

% -------------------------------------------------------------------------
if nargin<5
    fi_shift = zeros(size(amp));
end

amp = amp(:)';
freq = freq(:)';
fi_shift = fi_shift(:)';
%

dR_sin = [];
time_vec=prt*(0:round(time/prt));
for A_f_s_ii = [amp;freq;fi_shift] % For number of targets
    amp_i  = A_f_s_ii(1);
    freq_i = A_f_s_ii(2);
    fi_i   = A_f_s_ii(3);
    dR_sin = [ dR_sin,...
                amp_i.*sin( 2*pi*freq_i*time_vec + fi_i )' ];
end



end