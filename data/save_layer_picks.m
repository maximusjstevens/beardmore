%%%%

% Picking layers for beardmore: run the px4i script. Go through and pick
% surface gets [xXXXs,yXXXs]=ginput_zoom (make sure that figure you want is
% current), where XXX is the last 3 numbers in the radar data filename. 
% 

% save a bunch of things:
num='813';
name=strcat(num,'picks.mat');
var1=strcat('x',num,'b'); % distance of the bed
var2=strcat('x',num,'s'); % distance surface
var3=strcat('y',num,'b'); % height bed
var4=strcat('y',num,'s'); % height surface

save(name,var1,var2,var3,var4,'Height_s','X','Y','Latitude','Longitude','Distance','Distance_total','HPOS')

% Also, save the figure that the picks were on: XXXpicks.mat

% For picks, x is Distance on track, y is the vertical position

% use Distance_total, Height_s

%First, interpolate picks onto common grid (use grid from radar survey)
xx=eval(var1);
yy=eval(var3);



if xx(1)>xx(end)
   %put flip lr here 

b_max=max(xx);
b_min=min(xx);

dist_interp=HPOS(:,3);

ind = find(dist_interp>=b_min);
ind2 = find(dist_interp<=b_max);

ind_bed = intersect(ind,ind2);

HPOS_bed = HPOS(ind_bed,:);

