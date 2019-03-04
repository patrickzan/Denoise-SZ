function [x,d]=nt_sns3(x,thresh)
% [y,w]=nt_sns2(x,thresh) - sensor noise suppression, new version
%
%  y: denoised data 
%  w: 0 for parts that needed fixing, 1 elsewhere (time*chans)
%
%  x: data to denoise (time*chans or time*chans*trials)
%  thresh: threshold for Mahalanobis distance (default:1);
% 
nt_greetings;

if nargin<1; error; end
if nargin<2 || isempty(thresh); thresh=1; end

[nsample,nchan,~]=size(x);
x=nt_unfold(x);

mn=mean(x); % save means
x=nt_demean(x);
nn=sqrt(mean(x.^2)); % save norms
x=nt_normcol(x);

%{
For each channel, find sections for which it does not fit the 
subspace spanned by other sensors.  The calculation is repeated 
and the projection matrix is refined at each step.
%}
PCA_THRESH=10^-15;

d=zeros(size(x));


% c0: covariance of non-artifact part
x=nt_demean(x);
c0=nt_cov(x);

for iChan=1:nchan

    oChan=setdiff(1:nchan,iChan); 

    % PCA to avoid problems with rank-deficient data
    [topcs,eigenvalues]=nt_pcarot(c0(oChan,oChan)); % PCA
    idx=find(eigenvalues/max(eigenvalues) > PCA_THRESH); % discard weak dims
    topcs=topcs(:,idx);

    % regression matrix for this channel on all others
    b=c0(iChan,oChan)*topcs/(topcs'*c0(oChan,oChan)*topcs);         
    y(:,iChan)=(x(:,oChan)*topcs)*b'; % projection 
    dd=y(:,iChan)-x(:,iChan); % difference from projection

    d(:,iChan)=mahal(dd,dd)/thresh; % excentricity of each sample from distribution

end    
%figure(1); clf; plot(d); pause

%{
Find which channel is most excentric at each moment.
%}
w=max(d,[],2)<1;
[~,iWorstChan]=sort(d','descend'); 
iWorstChan=iWorstChan';
iWorstChan(find(w),:)=nan;

x=nt_demean(x,w);
c0=nt_cov(x(find(w),:)); % covariance of non-artifact part
for iChan=1:nchan
    
    idx=find(iWorstChan(:,1)==iChan);
    
    ww=zeros(size(w));
    ww(idx)=1;
    
    if numel(idx)>0
        c1=nt_cov(x(idx,:));
        todss=nt_dss0(c0,c1);
        z=nt_mmat(x,todss);
        fromdss=pinv(todss);
        x(idx,iChan)=nt_mmat(z(idx,2:end),fromdss(2:end,iChan));
    end
    

end


x=bsxfun(@times,x,nn);
x=bsxfun(@plus,x,mn);

x=nt_fold(x,nsample);
d=nt_fold(d,nsample);

