% f_gps.m
% twit, 2013
% routine adapted from f_press (trades, 1997) to correct the surface elevation using GPS data and PXI receiver
%
% This routine takes the filtdata array and shifts each record "vertically" a distance
%  based on the GPS. 
% Assume pretrigger data in filtdata have already been removed 
%  The shift of each record is scaled to ice thickness using the "Cice" value.  
%
% The GPS record is lowpass filtered first.

 Cice=(169e6)/2;

% smooth the GPS Record using a lowpass filter hard coded at 500m;
  inc=abs(mean(diff(X_s)));
%       if(inc<1)  % then X_s is in Km; we want meters
%	inc=inc*1000;
%      end

   Sample_Freq = 1/(inc);  	%inc is in units: meters/record

    Nyquist_Freq = 1/2*Sample_Freq;

    Lowpass_Freq = 1/500;  % smooth below 500m

    Corner_Freq = Lowpass_Freq/Nyquist_Freq;

% calculate filter polynomials for 4th order butterworth lowpass

      [b,a] = butter(4, Corner_Freq);

%  sweep through GPS record and filter it
%

       Height_s1 = filtfilt(b, a, Height_s);

  clear inc Sample_Freq Nyquist_Freq Lowpass_Freq Corner_Freq b a 

sz = size(filtdata);
   Vpos=Height_s1;
   Displacement = 1 + round( Vpos/Cice*1/(Sampling_int*1e-9));
   max_disp = abs(max(Displacement));

sz(1) = sz(1)+max_disp+1;

sz2=sz;

sz=size(filtdata);

% expand the filtdata array to the size necessary to accomodate the most shifted
%   array.
filtdata=[filtdata ; zeros(sz2(1)-sz(1),sz(2))];


% now figure out the section specific to this profile
clear Vpos Displacement

Vpos = Height_s1;
Displacement = 1 + round(Vpos/Cice*1/(Sampling_int*1e-9));


for i = 1:sz(2)
        filtdata( 1+max_disp-Displacement(i):max_disp-Displacement(i)+sz(1),i) = filtdata(1:sz(1),i);
        filtdata(1:1+max_disp-Displacement(i)-1, i ) = zeros(size(filtdata(1:1+max_disp-Displacement(i)-1, i )));  
end










sz = size(filtdata);

TWTtime_old=TWTtime;

for i=1:sz(1)
	TWTtime(i) = (i)*(Sampling_int)*1e-9;
end

