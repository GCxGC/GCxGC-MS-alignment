function [rtpix] = Time2Pix(rttime,Mod,Freq,isot)
% From rtpix a retention vector (in times), a point per line 
% (input in (s;min), output in (pixel;pixel) (with first value at (1;1))
% Mod is modulation rate, in [s], Freq is sampling frequency, in [Hz], isot
% =if you want to suppress some time at the beginning of the chromatogram,
% in [min].
% (by default, isot = 0 if not given)
if ~exist('isot','var')
    isot = 0;
end
rtpix(:,2)=round(rttime(:,2)*Freq);
rtpix(:,1)=round((rttime(:,1)-isot)*60/Mod);