%  THIS FILE READS DATA FROM CNT INTO MATLAB AND DOES 
% 1. DO A BASIC HIGHPASS FILTERING. 
% 2. REJECT BAD CHANNELS..SO THAT IT DOESN'T EFFECT RE-REREFERENCING.
% 3. AVERAGE RE-REFERENCING IF NO OF CHANNELS > 32  
%       ? should referencing include ocular channels ??
% 4. REJECT BAD EPOCHS AFTERS TSPCA.


function [fs,data_highpass, labels_electrodes, bad_channels,active_channels,events_trigger] = ...
    cnt_matlab_conversion(CNT_file_name, fs, CNT_file_directory, cc)

bits_data='int32';

[f, lab, dat, labels_electrodes, events_trigger, ~] = Read_Neuroscan_COGMO(CNT_file_name,bits_data,CNT_file_directory);


%% high pass filter at 0.8 hz before avg re-ref: IIR butter filter
if cc == 1
    Lcut = 0.4;
else
    Lcut = 0.8;
end
Wn = Lcut/(fs/2);                   % Normalized cutoff frequency at 0.8 Hz.      
[z,p,k] = butter(15,Wn,'high');  % Butterworth filter
[sos,G] = zp2sos(z,p,k);          % Convert to SOS form
%h = fvtool(sos);

data=dat';  % data is time X channels
data_highpass=zeros(size(data));

for i=1:size(data,2)
    data_highpass(:,i)=filtfilt(sos,G,data(:,i));   %this is data highpass at 1hz.
end

%% notch filter at 60 ,120 and 180.
Q=60; Wo = 60/(fs/2);  BW = Wo/Q;[b,a] = iirnotch(Wo,BW);
[sos, G]=tf2sos(b,a);

for i=1:size(data_highpass,2)
    data_highpass(:,i)=filtfilt(sos,G,data_highpass(:,i));   %60hz notch.
end

Q=60; Wo = 120/(fs/2);  BW = Wo/Q;[b,a] = iirnotch(Wo,BW);
[sos, G]=tf2sos(b,a);
for i=1:size(data_highpass,2)
    data_highpass(:,i)=filtfilt(sos,G,data_highpass(:,i));   %120hz notch.
end


Q=60; Wo = 180/(fs/2);  BW = Wo/Q;[b,a] = iirnotch(Wo,BW);
[sos, G]=tf2sos(b,a);
for i=1:size(data_highpass,2)
    data_highpass(:,i)=filtfilt(sos,G,data_highpass(:,i));   %180hz notch.
end


ocular_channels=[];
ocular_channels=[ocular_channels ; find(strcmp('HEO',cellstr(labels_electrodes)))];
ocular_channels=[ocular_channels ; find(strcmp('VEO',cellstr(labels_electrodes)))];
ocular_channels=[ocular_channels ; find(strcmp('HEOG',cellstr(labels_electrodes)))];
ocular_channels=[ocular_channels ; find(strcmp('VEOG',cellstr(labels_electrodes)))];
active_channels=setdiff([1:length(labels_electrodes)],ocular_channels);
%% reject bad channels.

%reject only among active channels.. we don't know how the ocular channel
%statistics are going to be.

bad_channels=[];

% variance criterion
channel_var=var(data_highpass(:,active_channels),[],1);
indices=find(channel_var > mean(channel_var)+3*std(channel_var));
bad_channels=union(indices,bad_channels);
%% change to average reference if no of electrodes > 32 otherwise leave it..
good_channels = setdiff(active_channels, bad_channels);

if size(dat,1) >32
    dat_mean = mean(data_highpass(:, good_channels),2);  %avg doesnot include ocualr channels.
    %plot_avg_concat_spectrum(dat_mean-mean(dat_mean),events_trigger,SSR_rate);
    data_highpass = data_highpass-repmat(dat_mean,1,size(data_highpass,2));
end

data_highpass = demean(data_highpass);

%kurtosis criterion
% there are few outlier data points which screwup the kurtosis and make its
% value veryhigh.. so lets not put kurtosis criterion for now...
%channel_kur=kurtosis(data_highpass(:,active_channels),[],1);


end
