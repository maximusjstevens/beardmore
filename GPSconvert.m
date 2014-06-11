
% input_id=input('Input filename (Year+Julian Date+Hour+Min+Sec): ', 's');
 id=strcat(input_id, 'GPS');

%id='2012335125735GPS';

fid=fopen(id);

GPSdata=textscan(fid, '%s %s %s %f %s %s %f %s %f %f %f %f %s');

fclose(fid);

clear ans fid

date=cell2mat(GPSdata{1});
time=cell2mat(GPSdata{2});

str=cell2mat(GPSdata{3});
digstr=str(:,1:2);
latD=str2num(digstr);

latM=GPSdata{4};
latM=latM/60;

lat=latD+latM;

latNS=cell2mat(GPSdata{5});

str=cell2mat(GPSdata{6});
digstr=str(:,1:3);
lonD=str2num(digstr);

lonM=GPSdata{7};
lonM=lonM/60;

lon=lonD+lonM;

lonEW=cell2mat(GPSdata{8});

Z=GPSdata{9};
latERR=GPSdata{10};
lonERR=GPSdata{11};
ZERR=GPSdata{12};

latsign=zeros(size(lat));
latPOS=find(latNS=='N');
latsign(latPOS)=1;
latNEG=find(latNS=='S');
latsign(latNEG)=-1;

lonsign=zeros(size(lon));
lonPOS=find(lonEW=='E');
lonsign(lonPOS)=1;
lonNEG=find(lonEW=='W');
lonsign(lonNEG)=-1;

coor=[lonsign.*lon latsign.*lat];
coor=coor';

savename=sprintf('coor%s.txt',id);
fileID=fopen(savename,'w');
fprintf(fileID,'%s\t %s\n','Lon','Lat');
fprintf(fileID,'%.6f \t %.6f \n',coor);
% save(savename,'coor','-ascii','-tabs')



