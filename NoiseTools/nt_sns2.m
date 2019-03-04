function [x,w]=nt_sns2(x,thresh)
% [y,w]=nt_sns2(x,thresh) - sensor noise suppression, new version
%
%  y: denoised data 
%  w: 0 for parts that needed fixing, 1 elsewhere (time*chans)
%
%  x: data to denoise (time*chans or time*chans*trials)
%  thresh: threshold for Mahalanobis distance (default:1);
% 
nt_greetings;

PCA_THRESH=10^-15;

if nargin<1; error; end
if nargin<2 || isempty(thresh); thresh=1; end

[nsample,nchan,~]=size(x);
x=nt_unfold(x);

mn=mean(x); % save means
x=nt_demean(x);
nn=sqrt(mean(x.^2)); % save norm
x=nt_normcol(x);

%{
For each channel, find sections for which it does not fit the 
subspace spanned by other sensors.  The calculation is repeated 
and the projection matrix is refined at each step.
%}

w=ones(size(x,1),1);

NITER=2; % iterations to refine c0
for k=1:NITER
    
    % c0: covariance of non-artifact part
    x=nt_demean(x,w);
    c0=nt_cov(x,[],w);
    
    for iChan=1:nchan

        % regression matrix for this channel on all others
        oChan=setdiff(1:nchan,iChan); 
        
        [topcs,eigenvalues]=nt_pcarot(c0(oChan,oChan)); % PCA
        idx=find(eigenvalues/max(eigenvalues) > PCA_THRESH); % discard weak dims
        topcs=topcs(:,idx);
        b=c0(iChan,oChan)*topcs/(topcs'*c0(oChan,oChan)*topcs); % matrix to project on other channels        
        y(:,iChan)=(x(:,oChan)*topcs)*b'; % projection 
        dd=y(:,iChan)-x(:,iChan); % difference from projection
        %plot([y(:,iChan),x(:,iChan)]); pause
        d=mahal(dd,dd)/thresh; % excentricity of each sample and each channel 
        w=min(w,(d<1));

    end    
    disp(mean(w))
end

x=nt_demean(x,w);
c0=nt_cov(x,[],w);


%{
For each channel, find the part for which it dominates the artifact.
Replace that part based on projection on the other channels. 
%}

%plot(w); pause

iBad=find(~w);
xx=x(iBad,:);
z=nt_pca(xx);
xx=abs(nt_normcol(xx));

%figure(1); clf; plot(xx); pause

%xx=nt_smooth(xx,20,[],1);

for iChan=1:nchan
    
    iBad2=find(xx(:,iChan)>max(xx(:,setdiff(1:nchan,iChan)),[],2)); % this channel dominates others
    oChan=setdiff(1:nchan,iChan);
    [topcs,eigenvalues]=nt_pcarot(c0(oChan,oChan)); % PCA        
    idx=find(eigenvalues/max(eigenvalues) > PCA_THRESH); % discard weak dims
    topcs=topcs(:,idx);
    b=c0(iChan,oChan)*topcs/(topcs'*c0(oChan,oChan)*topcs); % matrix to project on other channels        
    x(iBad(iBad2),iChan)=(x(iBad(iBad2),oChan)*topcs)*b'; % projection 
end

%{
To do:
Record all the DC shifts introduced when x mean is removed, so as to
restore accurately.
%}


x=nt_demean(x);
x=bsxfun(@times,x,nn);
x=bsxfun(@plus,x,mn);

x=nt_fold(x,nsample);
w=nt_fold(w,nsample);



