function [x,y]=lltops(lat,lon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% turns lat,lon pairs into polar stereographic pairs
% uses E=.082271850
% R=6378.2650 km
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  


%constants:
  PI  =	4*atan(1);
  DR  = PI/180;
  CDR = 180/PI;
  %RE  = 6378206.4;
  RE=6378137.0;  % for WGS84
  %E2  = 0.006768658;
  %E   = 0.082271850;
  E=0.08199188997903; % for WGS84;
  E2=0.00672267002233; % for WGS84;
  SLAT= 71.0;
  t   = tan(PI/4.0 - SLAT/(2.0*CDR)) ./ ( ((1.0 - E*sin(DR*SLAT))/ (1.0 + E*sin(DR*SLAT)))^(E/2.0) );
  cm  =	cos(DR*SLAT)/sqrt(1.0 - E2*(sin(DR*SLAT)*sin(DR*SLAT)));

  sn=1.0;
  sn=sign(lat);
        lat = lat .* sn;
        t1 = tan ( PI / 4.0 - lat / ( 2.0 * CDR) )./( ( ( 1.0 - E * sin( DR * lat ) ) ./ ( 1.0 + E * sin( DR * lat ) ) ) .^ ( E / 2.0 ));
        rho = RE * cm * t1 / t;
        x =  real(-rho.*sn.*sin(DR*lon));
        y = real(rho.* cos(DR*lon));





