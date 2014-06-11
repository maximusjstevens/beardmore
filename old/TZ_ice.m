function[z_out] = TZ_ice(t_in);  
% t_in is 2-way travel time as seen as Rx, with time zeroed on airwave arrival 
% z_out is the corresponding relfection depth   
c_ice = 1.685e8; 
c = 3e8;  
% distance from Tx-sled to Rx-sled  
SEP = 126;  
% antennas are parallel and air wave distance is shortened by DEL_SEP 
DEL_SEP = 0;  
% ray geometry 
trav_dist = c_ice * (t_in + (SEP-DEL_SEP)/c);  
z_out = sqrt( max( (trav_dist/2).^2 - (SEP/2)^2 , 0 ) );