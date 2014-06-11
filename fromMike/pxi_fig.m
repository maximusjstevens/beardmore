% ReadPXIdata.m % This script is to read radar data collected by PXI system.  % The data format was defined at August 23, 2005.  % GPS always sends time info at UTC.  If GPS worked, both GPS time stamp % and radar time stamp are given by the GPS.  However, if GPS didn't work,
% radar time stamp corresponds to the PXI system time, which may be UTC or
% local time.  Anyway, the radar operator should correct PXI time and host
% computer time with a handheld GPS.  It prevents a lot of confusion in the
% data processing.  This script assumes that PXI system time is set to the
% local time.

% +++++ Revision record ++++++
% === November 29, 2005
% output filename was revised.  Add if-end for GPS processing, for the case
% if the GPS is not available.  McMurdo time zone was revised to -11 from
% -12 because we go down in the summer only.  This revised one will be
% used Inland WAIS divide traverse.
% === December 9, 2005
% gps location map is now axis "equal".  Data summary shows GPS data
% aquatision successful rate, rather than collapsed rate.  UTC to local
% time conversion was revised so that now the correct local time can be
% derived.
% === December 10, 2005
% Bug in TWTtime is fixed.
% === December 12, 2005
% Bug in the file saving is fixed.  Now, multiple-line notes can be read.
% Numeric value of Timezone_flag is now replaced by a string so that its
% mean is more apparent.  This replacement is done at the last part of this
% program.  Fid is now cleared at the end of the script.
% === January 17, 2006
% Filename error message is revised.  Graph colors for hour/min/sec was
% revised.
% === January 17, 2007
% New figure about data collection timing was added.

dir 20*info
Filename = input('Input filename (Year+Julian Date+Hour+Min+Sec): ','s');
Filename_info = strcat(Filename,'info');
Filename_gps = strcat(Filename,'gps');
Filename_num = strcat(Filename,'num');
Filename_notes = strcat(Filename,'notes');
Filename_ch0 = strcat(Filename,'ch0');
Filename_ch1 = strcat(Filename,'ch1');



%===== Check all files are now copied to the target directory
%Filename_array = char(Filename_info, Filename_gps, Filename_num, Filename_notes, ...
%    Filename_ch0, Filename_ch1);
%It's inserting a space for some reason. remove spaces
%Filename_array=Filename_array(:,(1:size(Filename_array,2)-1))
%
%for I = 1:6
%    Fid = fopen(Filename_array(I,:), 'r');
%    if Fid == -1
%        switch I
%            case {1,3,4,5}
%                Filename_array(I,:)
%                error( strcat(Filename_array(I,:), ': unavailable.  File name might be typed incorrectly.') );
%            case {2}
%                disp('GPS data are not avaialble');
%            otherwise
%                disp('ch1 data are not available.  Only ch0 data will be read');
%        end
%    else
%        fclose(Fid);
%    end
%end
%clear Filename_array;
    
%===== Read notes file.
%Fid_notes = fopen(Filename_notes, 'rt')
%if feof(Fid_notes)
%    Notes(1, :) = 'no notes';
%else
%    I = 1;
%    Length_notes = 0;
%    while feof(Fid_notes) ~= 1 % count lines in the notes.
%       fgetl(Fid_notes);
%        I = I+1;
%    end
%    Length_notes = I-1;
%    fseek(Fid_notes, 0, 'bof');
%    Notes = num2str( fgetl(Fid_notes) );
%    for I = 2:Length_notes
%        Temp = num2str( fgetl(Fid_notes) );
%        Notes = char(Notes, Temp);
%        clear Temp;
%    end
%    fclose(Fid_notes);
%    clear I;
%end

% ===== Display notes and make sure this is the file that the user wants to
% read into
%disp('===== Notes =====');
%disp(Notes);
%disp('=================');
%Flag_notes = input('Do you want to read this file? If yes, input 1:');
%if Flag_notes ~= 1
%    return;
%end

%Timezone = input('Time difference from the UTC (McMurdo: -11, Alaska summer time: -8)? ');
Timezone = -13;
%===== Read info file.
[GPS_onoff, Ref_point, Sampling_int, Datalength, Stacking...
    Ch0, Range_0, Offset_0, Ch1, Range_1, Offset_1, ...
    Impedance, Trig_source, Trig_coupling] ...
    = textread(Filename_info, '%s Refpoint%d SampInt%f Datalength%d Stacking%d Ch%d %f %f Ch%d %f %f Impedance%d %s %s', 1);


% Info file sample format
% gpsON Refpoint5 SampInt20 Datalength2500 Stacking1000 Ch0 2 -1 Ch1 0 0 Impedance50 EXT ac
%         %          nsec                               range offset     Impedance, trigger source, coupling
clear Ch0 Ch1;

%===== Check file size consistency
Num = textread(Filename_num, '%d', 1);
if Num <= 0
    error('Num == 0');
    return
end
%===== GIVE NUM HERE, ONLY FOR TESTING PURPOSE =====%
% Num = 100;
% This script assumes that Num is the most reliable number that shows the
% number of data sets.  When the Labview program starts, 10 of zero is
% padded in the file.  Num counter is incremented every time when the new
% data set was collected, but the file is not updated. When the all files
% are written successfuly or the measurement is terminated due to errors,
% the Num ib the file is updated.


%===== Read waveform data from each channel
if Range_0 ~= 0 % 0 stands for no data collection from Ch1
    Fid_ch0 = fopen(Filename_ch0, 'r');
    Temp_ch0 = fread(Fid_ch0, [Datalength+3 Num], 'double', 'b'); % ieee-big endian option
    fclose(Fid_ch0);
end
if Range_1 ~= 0 % 0 stands for no data collection from Ch1
    Fid_ch1 = fopen(Filename_ch1, 'r');
    Temp_ch1 = fread(Fid_ch1, [Datalength+3 Num], 'double', 'b'); % ieee-big endian option
%    Temp_ch1 = fread(Fid_ch1, [Datalength Num], 'double', 'b'); % ieee-big endian option
    fclose(Fid_ch1);
end

%===== This comment part is an old script to check the size of data file
%and Num.  However, the policy was changed so that we now believe Num.
% %===== Check file size consistency
% Num = textread(Filename_num, '%d', 1);
% M0 =0; M1 = 0; N0 = 0; N1 = 0;
% if Range_0 ~= 0
%     [M0, N0] = size(Temp_ch0);
% end
% if Range_1 ~= 0
%     [M1, N1] = size(Temp_ch1);
%     if N0 ~= N1
%         error('Data size for ch0 and ch1 is NOT consistent.');
%     end
% end
% 
% if Num == 0
%     Num = max(N0, N1); % 'max' works, if only one channel data were collected.
%     disp('Num file is collapsed.  Data file size was used to define the number of waveform.');
% elseif Num ~= max(N0, N1)
%     error('Num and data size is NOT consistent.');
%     return;    
% end
% clear M0 M1 N0 N1;

%===== Read GPS data and then convert to the local time
if strcmp(GPS_onoff, 'gpsON')
    [MonthU, DateU, YearU, HourU, Minute, Second, ...
        Latitude_deg Latitude_min, SN, ...
        Longitude_deg, Longitude_min, EW, ...
        Height, HPE, VPE, EPE GPS_status] = ...
        textread(Filename_gps, '%d/%d/%d %d:%d:%d %ddeg %f %s %ddeg %f %s %f %f %f %f %s', Num);
    % GPS sample format
    % 08/24/2005 16:50:20 47deg 39.2456 N 122deg 18.5544 W -18.4 7.7 8.4 11.4 Valid<CR><LF>
    % MO/DATE/YEAR  TIME Latitude Longitude Height HPE VPE EPE GPS_status
    % This height is geoid height, not elevation above sea level.
    % If serial port communication is incorrect, the string will be
    % replaced by
    % 99/99/2099 99:99:99 0deg 0 N 180deg 0 W -999 -1 -1 -1 Collapse<CR><LF>
    
    for I = 1:Num % If GPS data is collapsed, replace the time stamp.
        if MonthU ~= 99
            Closest = I;
            I = I+1;
        else
            MonthU(I) = MonthU(Closest);
            DateU(I) = DateU(Closest);
            YearU(I) = YearU(Closest);
            HourU(I) = HourU(Closest);
            Minute(I) = Minute(Closest);
            Second(I) = Second(Closest);
            I = I+1;
        end
    end
    
    % UTC -> Local time conversion
    HourL = HourU + Timezone;
    Timezone_flag = 1; % Flag on for message when done.
    if HourL > 24 & HourL < 37
        DateL = DateU - 1;
        HourL = HourL - 24;
    elseif HourL < 0
        DateL = DateU + 1;
        HourL = 24 + HourL;
    else
        DateL = DateU;
        Timezone_flag = 0; % Flag off.
    end
    YearL = YearU; MonthL = MonthU;
    
    % Latitude/Longitude conversion
    Latitude = Latitude_deg + Latitude_min./60;
    if strcmp(SN, 'S')
        Latitude = -1 .* Latitude;
    end
    Longitude = Longitude_deg + Longitude_min./60;
    if strcmp(EW, 'W')
        Longitude = -1 .* Longitude;
    end
    clear Latitude_deg Latitude_min Longitude_deg Longitude_min;
end

%===== Divide ch0/1 matrix to radar waveform and time stamp
% If GPS is not availanle, HourR shows the local time, not UTC.
%if Range_0 ~= 0 % Grab radar time stamp from either Ch0 or Ch1
%    HourR = Temp_ch0(1, 1:Num)';
%    MinuteR = Temp_ch0(2, 1:Num)';
%    SecondR = Temp_ch0(3, 1:Num)';
%elseif Range_1 ~= 0
    HourR = Temp_ch1(1, 1:Num)';
    MinuteR = Temp_ch1(2, 1:Num)';
    SecondR = Temp_ch1(3, 1:Num)';
%end
if Range_0 ~= 0
    size(Temp_ch0)
    Num
    Waveform_ch0 = Temp_ch0(4:Datalength+3, 1:Num);
    clear Temp_ch0;
end
if Range_1 ~= 0
    Waveform_ch1 = Temp_ch1(4:Datalength+3, 1:Num);
    clear Temp_ch1;
end

if strcmp(GPS_onoff, 'gpsON')
    %===== Generate polar stereographic coordinate
    % No check will be made whether polar stereographic projection is appropriate or not.
    % UTM coordinate will be added to this routine in the future.
    [X, Y] =ll2ps(Latitude, Longitude);
    %Distance = zeros(Num-1, 1);
    %Distance_total = zeros(Num-1, 1);
    %    for I = 0 : Num-1
    %    Distance(I) = sqrt( (X(I+1)-X(I)).^2 + (Y(I+1)-Y(I)).^2 );
    %    if I == 0
    %        Distance_total(I) = Distance(I);
    %    else
    %        Distance_total(I) = Distance_total(I-1) + Distance(I);
    %    end
    %end
    %clear I;

    Distance = zeros(Num, 1);
    Distance_total = zeros(Num, 1);
       for I = 2 : Num
        Distance(I) = sqrt( (X(I)-X(I-1)).^2 + (Y(I)-Y(I-1)).^2 );
            if Distance(I)<1;
        Distance(I)=0;
            end
        Distance_total(I) = Distance_total(I-1) + Distance(I);
            if Distance_total(I-1)==0;
                Istart=I-1;
            end
%            if Distance_total(I)>Distance_total(I-1);
                Iend=I+1;
%            end
       end         
    clear I;
    

% =============Before filtering
    %1 =Remove pretrigger rows
%Waveform_ch0(1:Ref_point*Datalength/100, :)=[];
Waveform_ch1(1:Ref_point*Datalength/100, :)=[];

    %2 =Remove cols when system stationary at start and end of profiles
%Waveform_ch0(: ,Iend:Num)= [];
%Waveform_ch0(: ,1:Istart-1)= [];
Waveform_ch1(:,Iend:Num) = [];
Waveform_ch1(: , 1:Istart-1) = [];

%demean ch1 waveforms along a transect
   [m n] = size(Waveform_ch1);
    mean_trace = mean(Waveform_ch1, 2); % mean trace along transect
    Waveform_ch1 = Waveform_ch1 - mean_trace(:, ones(1, n));

    %=====filtering ch1, a la Ben
  [b,a]=butter(3, [.3 4]/25); 
  filtdata=filtfilt(b, a, detrend(Waveform_ch1));
 % [b,a]=butter(3, [.5 4]/25);
 %filtdata2=filtfilt(b, a, detrend(Waveform_ch1));
 %[b,a]=butter(3, [.4 4]/25); 
 % [b,a]=butter(3, [5 12]/25);
  %  fdata=filtfilt(b, a, detrend(convstack(Waveform_ch1, 4)));  % uncomment
   % to stack by 4
 % filtdata3=filtfilt(b, a, detrend(Waveform_ch1));  
 % filtdata=filtfilt(b, a, (Waveform_ch1));
 %clear Waveform_ch1


   
%demean ch0 waveforms along a transect
%    [m n] = size(Waveform_ch0);
%    mean_trace_0 = mean(Waveform_ch0, 2); % mean trace along transect
%    Waveform_ch0 = Waveform_ch0 - mean_trace_0(:, ones(1, n));

    %=====filtering ch0, a la Ben
  %[b,a]=butter(3, [.5 8]/25); 
  %[b,a]=butter(3, [4 12]/25);
%    [b,a]=butter(3, [5 12]/25);
    %fdata=filtfilt(b, a, detrend(convstack(Waveform_ch1, 4)));  % uncomment
   % to stack by 4
%   filtdata_0=filtfilt(b, a, detrend(Waveform_ch0));   
%==== Generate time and depth scale
% For low frequency ground-based measurements (no firn). Pretrigger already
% removed from waveforms
TWTtime = 0:Datalength-(Ref_point*Datalength/100)-1;
TWTtime = (Sampling_int .* 1e-9) .* TWTtime; %Sampling_int is given at nsec.
%TWTtime(1:Ref_point*Datalength/100, :)=[];
t_in = TWTtime;
z_out = TZ_ice(t_in);
clear TWTtime;
clear t_in;

%==== Adjust depth matrix
%z_out(1:11)=[];
%filtdata(1:11, :)=[];
%==== Adjust horizontal distance
Distance_total(Iend:Num)= [];
Distance_total(1:Istart-1)= [];

%===== Generate voltage scale
% Correct reading above gives the real voltage scale.  It does not require
% any conversion using voltage full-scale and shift.

% Generate screened GPS data
    X_s = X; Y_s = Y; Height_s = Height; HPE_s = HPE; VPE_s = VPE; EPE_s = EPE;

    X_s( find( strcmp(GPS_status, 'Collapse') ) ) = NaN;
    Y_s( find( strcmp(GPS_status, 'Collapse') ) ) = NaN;
    Height_s( find( strcmp(GPS_status, 'Collapse') ) ) = NaN;
    HPE_s( find( strcmp(GPS_status, 'Collapse') ) ) = NaN;
    VPE_s( find( strcmp(GPS_status, 'Collapse') ) ) = NaN;
    EPE_s( find( strcmp(GPS_status, 'Collapse') ) ) = NaN;

    %==== Adjust screened height matrix
Height_s(Iend:Num)= [];
Height_s(1:Istart-1)= [];

    % Remove NaN so that derive static values
    Temp = Height_s; Temp( find(isnan(Height_s)) ) = [];
    Height_mean = mean(Temp);
    Height_max = max(Temp);
    Height_min = min(Temp);
    clear Temp;
    Temp = HPE_s; Temp( find(isnan(HPE_s)) ) = [];
    HPE_mean = mean(Temp);
    HPE_dev = std(Temp);
    clear Temp
    Temp = VPE_s; Temp( find(isnan(VPE_s)) ) = [];
    VPE_mean = mean(Temp);
    VPE_dev = std(Temp);
    clear Temp;
    Temp = EPE_s; Temp( find(isnan(EPE_s)) ) = [];
    EPE_mean = mean(Temp);
    EPE_dev = std(Temp);
    clear Temp;
end

%===== Draw graphs to check data quality

%figure; %radargram  Fig. 1
%if (Range_0 ~= 0) & (Range_1 ~= 0)
%    subplot(2,1,1);
%    imagesc(Waveform_ch0); grid on;
%    colorbar;
%    title('Ch0 Radar Data');
%    subplot(2,1,2);
%    imagesc(Waveform_ch1); grid on;
%    title('Ch1 Radar Data');
%    colorbar;
%elseif Range_0 ~= 0
%    imagesc(Waveform_ch0); grid on; title('Ch0 Radar Data');
%    colorbar;
%else
%    imagesc(Waveform_ch1); grid on; title('Ch1 Radar Data');
%    colorbar;
%end

%if strcmp(GPS_onoff, 'gpsON') % difference between GPS time stamp and radar time stamp
%    figure; % time stamp check, Fig. 2
%    subplot(2,2,1);
%    title('Local time stamp plot');
%    plot(DateL, 'y'); grid on; legend('Date');
%    subplot(2,2,2);
%    plot(HourL, 'b'); grid on;  legend('Hour');
%    subplot(2,2,3);
%    plot(Minute, 'r'); grid on;  legend('Minute');
%    subplot(2,2,4);
%    plot(Second, 'r'); grid on;  legend('Second');
%
%    figure; % Fig. 3
%    plot(MinuteR - Minute, 'r'); hold on;
%    plot(SecondR - Second, 'b');
%    grid on;
%    legend('Minute difference', 'Second difference');
%    title('Radar time stamp - GPS time stamp (it must be zero)');
%
%    Flag = input('Do you want to map the profile location? If yes, input 0: ');
%    if Flag == 0
%%%        figure; % Fig. 4
%%%        subplot(2,1,1);
        %plot(X./1000, Y./1000, '-bp'); hold on;
%%%        plot(X_s./1000, Y_s./1000, '-rp');
        T = 0:pi./100:2.*pi;
        %plot( (cos(T).*HPE_mean./2+X(1))./1000, (sin(T).*HPE_mean./2+Y(1))./1000, 'r:');
        %plot( (cos(T).*(HPE_mean+HPE_dev)./2+X(1))./1000, (sin(T).*(HPE_mean+HPE_dev)./2+Y(1))./1000, 'g:');
        clear T;
%11        axis equal;
%11        grid on;
        %legend('raw', 'screened', 'mean error', 'mean+dev error');
%11        title(strcat('Map profile: mean error = ', num2str(mean(HPE), '%3.1f'), ' m.'));
%11        xlabel('(km)');
%11        ylabel('(km)');

 %1       figure; % Fig. 5
 %11       subplot(2,1,2);
 %11       plot(Distance_total/1000, Height_s);
 %11       grid on;
 %11       xlabel('(km)');
%%        ylabel('Elev above WGS84 (m)');
        %hold on;
        %plot(Height_s + VPE_mean, 'b:');
        %plot(Height_s - VPE_mean, 'b:');
        %plot(Height_s + (VPE_mean + VPE_dev), 'g:');
        %plot(Height_s - (VPE_mean + VPE_dev), 'g:');
        %title('Elevation above WGS84');
    %end
    %clear Flag

%else

figure; imagesc (Distance_total/1000, z_out, filtdata, [-.4/100 .4/100]); grid on;  
colormap (bone(256)); axis([0 max(Distance_total/1000) 0 max(z_out)]);

figure; grid on; pcolor (Distance_total/1000, -z_out, filtdata); grid on;  
colormap (bone); caxis([-.5/100 .5/100]); shading flat;

figure; imagesc (Distance_total/1000, z_out, filtdata, [-.08/100 .08/100]); grid on;  
colormap (bone(256)); axis([0 max(Distance_total/1000) 0 max(z_out)]);

figure; imagesc (Distance_total/1000, z_out, filtdata, [-3/100 3/100]); grid on;  
colormap (bone); axis([0 max(Distance_total/1000) 0 3000]);

figure; imagesc (Distance_total/1000, z_out, Waveform_ch1, [-3/100 3/100]); grid on;
colormap (bone); axis([0 max(Distance_total/1000) 0 max(z_out)]);

Convert_time = fix(clock);
clear Filename_* ans

disp('%%% Radar data were converted to a matlab file successfuly %%%');
disp('===== SUMMERY =====:');
disp( strcat('Number of radar data: ', num2str(Num)) );
%disp( strcat('Equivalent maximum ice thicknbess: ', num2str(Depth(Datalength), '%5.0f'), ' (m)' ));
disp( strcat('Stacking number: ', num2str(Stacking, '%5.0f')) );
if strcmp(GPS_onoff, 'gpsON')
    disp( strcat('GPS status: ON') );
    if Timezone_flag == 1
        disp('-----');
        disp('UTC -> Local time conversion causes date change.');
        disp('Also, script simply copied YearU to yearL, and MonthU to MonthL');
        disp('Check time stamp carefully, and make a manual change if necessary'); 
        disp(strcat('Dates ranges from ', int2str(DateL(1)), ' to ', int2str(DateL(Num)), '.'));
        disp('-----');
        Timezone_flag = 'UTC -> Local time conversion causes date change.';
    else
        disp('UTC -> Local time conversion does not change the date.');
        Timezone_flag = 'UTC -> Local time conversion DOES NOT cause date change.';
    end
%    disp( strcat('Total distance of radar profile: ', num2str(Distance_total(Num-1)./1000, '%1.1f'), ' (km)') );
%   disp( strcat('Distance interval between radar data: ', num2str(Distance_total(Num-1)./(Num-1), '%1.2f'), '(m)') );
    disp( strcat('Height difference: ', num2str(Height_max - Height_min, '%4.0f'), ' (m) with mean error of', num2str(VPE_mean, '%4.0f'), ' (m)') );
    disp( strcat('GPS data successful rate: ', num2str(length(isnan(X))./Num.*100, '%2.1f'), ' (%)') );
else
    disp( strcat('GPS status: OFF') );
    disp('GPS measurements are not available and Radar time stamp is given by the system clock.');
    disp('Check PXI system time clock and time zone (NO UTM time stamp is available).');
end
disp( strcat('This conversion was done at: (yyyy mm dd hh mm ss) ', num2str(Convert_time), ' (local time)') );

Flag1 = input( strcat('Do you want to save this file? If yes, input 0: ') );
if Flag1 == 0
    Fid = fopen(strcat(Filename, '.mat'), 'r');
    if Fid ~= -1
        fclose(Fid);
        Flag2 = input(strcat('Matlab file with the same name =', num2str(Filename), '.mat= is exist.\n',...
            'Save anyway with this file name, input 0.\n',...
            'Save the file with the time samp of file conversion, input 1.\n', ...
            'Cancel this save, input 2\n ??'));
            Flag2 = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            %%%%%%%%%%%%
        switch Flag2
            case 0
                save(Filename);
            case 1
                Filename = strcat(num2str(Convert_time(1)), ... % year
                    num2str(datenum(date) - datenum(YearL(1)-1, 12, 31)), ... % Julian day
                    num2str(Convert_time(4)), num2str(Convert_time(5)), ... % hour and minute
                    'converted_temp');
                save(Filename);
                disp('The data were saved into the file: YearJuliandayHourMinute_converted_temp');
            case 2
                return;
        end
        clear Flag2;
    else
        save(Filename);
    end
end

clear Flag1 Fid*;
