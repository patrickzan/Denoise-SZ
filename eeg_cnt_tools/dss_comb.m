function [to_dss, from_dss,ratio,pwr]=dss_comb(input,SSR_rate,fs)


%input is time X trials


%Mshifts is necessary to find the input data
%A bandpass filter is gonna to process the data
%The passband is between freqL and freqH

%In this version of DoDSS, both autocorrelation matrices are calculated as
%the mean of the autocorrelation matrices in different experiment
%conditions
%created by Nai 03/09
%calculate the covariance matrix block by block, modified by Nai 06/23/09
%based on the trial rejection code to select channels, modified by Nai 07/03/09

% S------filter design
% fs=100;SSR_rate=2.5;
%  q_factor=40;bw=(SSR_rate/(fs/2))/q_factor;
%  [b,a] = iircomb(fs/SSR_rate,bw,-6,'peak');  
%  b(end)=-1*b(end); a(end)=-1*a(end);  %to shift the peak positions
%   fvtool(b,a);
% 
%  [b1,a1]=butter(3,SSR_rate*4/(fs/2)); %3rd order lowpass filter.
%  fvtool(b,a,b1,a1);
%  
%  
%  Wn = 2/(fs/2);                   % Normalozed cutoff frequency        
% [z,p,k] = butter(15,Wn,'high');  % Butterworth filter
% [sos,G] = zp2sos(z,p,k);          % Convert to SOS form
% h = fvtool(sos);


% E------filter design



%prepare data
%and also calculate autocorrelation matrices
% cmat1 is the sphering autocorrelation matrix
% cmat2 is the biased autocorrelation matrix

%bad data rejection
sumch=squeeze(sum(abs(input(:,:,:)),2));
threshold=3*mean(mean(sumch,2));   %3 times the mean 
for ind=1:size(sumch,2)
    input(sumch(:,ind)>threshold,:,ind)=0;
end


%filtering data before giving to dss.
%not a good idea to demean per-trial basis.
% for trial=1:size(input,3)
%     input(:,:,trial)=demean(squeeze(input(:,:,trial)));
% end

%% fft method

% 
% for trial=1:1:size(input,3)
%     for channel=1:1:size(input,2)
%         input_filtered(:,channel,trial)=get_dss_bias(squeeze(input(:,channel,trial)),SSR_rate,fs);
%     end
% end
% 
% clean=input_filtered;
% % clean=clean(:,setdiff([1:157],MainSqd.badch),:);
% %large value rejection
% %     clean(abs(clean)>200)=0;
% 
% %cmat1
% inducedclean=unfold(input);
% cmat1(:,:)=inducedclean'*inducedclean;
% 
% %cmat 2
% evokedclean=mean(clean,3);
% cmat2(:,:)=evokedclean'*evokedclean*size(clean,3)^2;

%% comb filter method.

% for trial=1:1:size(input,3)
%     for channel=1:1:size(input,2)
%         input(:,channel,trial)=filtfilt(b,a,squeeze(input(:,channel,trial)));
%         input(:,channel,trial)=filtfilt(b1,a1,squeeze(input(:,channel,trial)));
%         input(:,channel,trial)=filtfilt(sos,G,squeeze(input(:,channel,trial))); %remove the first teeth of combfilter.
% 
%     end
% end

%% second order resonator at SSR rate with constant Q

 Wo = SSR_rate/(fs/2);  q_factor=40;BW = Wo/q_factor;
[b,a] = iirpeak(Wo,BW);
% fvtool(b,a, 'fs', fs);
unfolded_input=unfold(input);

for channel=1:size(unfolded_input,2)
    filtered_input(:,channel)=filtfilt(b,a,unfolded_input(:,channel));
end

%cmat1
inducedclean=unfolded_input;
cmat1(:,:)=inducedclean'*inducedclean;

% %cmat 2 for non-averaged
% filteredclean=filtered_input;
% cmat2(:,:)=filteredclean'*filteredclean;

%cmat 2 for averaged
avg_filtered=mean(fold(filtered_input,size(input,1)),3);
cmat2(:,:)=avg_filtered'*avg_filtered;


%% apply DSS using the average over epochs as a bias function
cmat2=mean(cmat2,3);cmat1=mean(cmat1,3);
keep2=10.^-13;
keep1=[];
% [todss,fromdss,ratio,pwr]=dss0(cmat1,cmat2,keep1,keep2);
[todss,fromdss,ratio,pwr]=dss0(cmat1,cmat2,keep1,keep2);
to_dss=todss;
from_dss=fromdss;
end




