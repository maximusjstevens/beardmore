function [lat,lon]=ps2ll(x,y);

if nargin==1 &  ~isreal(x)
   y=imag(x); x=real(x);
end

  PI	=	3.1415926535;
  DR    = 0.01745329251994;
  CDR   = 57.29577951;
  RE    = 6378206.4;
  E2    = 0.006768658;
  E     = 0.082271850;
  SLAT  = 71.0;
  t 	= 	tan(PI/4.0 - SLAT/(2.0*CDR)) / ( ((1.0 - E*sin(DR*SLAT))/(1.0 + E*sin(DR*SLAT)))^(E/2.0) );
 cm 	= 	cos(DR*SLAT)/sqrt(1.0 - E2*(sin(DR*SLAT)*sin(DR*SLAT)));

  sn=-1.0;
  xlam=180.0;

  rho = sqrt( x .* x + y .* y );
  chi = PI/2.0 - 2.0 * atan2( t * rho, RE * cm );
  lat = chi +...
     ((E2/2.0) + (5.0*E2*E2/24.0) + (E2*E2*E2/12.0))*sin(2.0*chi) + ...
     ((7.0*E2*E2/48.0) + (29.0*E2*E2*E2/240.0))*sin(4.0*chi) + ...
             (7.0*E2*E2*E2/120.0)*sin(6.0*chi);
  lat = lat*(sn * CDR);
  lon = xlam + sn * CDR * atan2(x, -y);
  lon((lon < 0.0))= lon(lon<0)+360.0;
  lon(lon > 360.0 )=lon(lon>360.0) - 360.0; 

