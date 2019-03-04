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
if nargin<1; error; end
if nargin<2 || isempty(thresh); thresh=1; end
PCA_THRESH=10^-15;
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

w=ones(size(x));

NITER=3; % iterations to refine c0
for k=1:NITER
    
    % c0: covariance of non-artifact part
    ww=min(w,[],2);
    x=nt_demean(x,ww);
    [c0,tw]=nt_cov2(x,ww);
    c0=c0./tw; c0(isnan(c0))=0;
    
    for iChan=1:nchan
        
        % regress this channel on all others
        oChan=setdiff(1:nchan,iChan); 
        
        [topcs,eigenvalues]=nt_pcarot(c0(oChan,oChan)); % PCA
        idx=find(eigenvalues/max(eigenvalues) > PCA_THRESH); % discard weak dims
        topcs=topcs(:,idx);
        b=c0(iChan,oChan)*topcs/(topcs'*c0(oChan,oChan)*topcs); % matrix to project on other channels
        y(:,iChan)=(x(:,oChan)*topcs)*b'; % projection 
        
        % fix parts where other channels are bad
        for iChan2=oChan
            oChan2=setdiff(oChan,iChan2);
            iBad=find(w(:,iChan2)==0);
            if ~isempty(iBad)
                [topcs,eigenvalues]=nt_pcarot(c0(oChan2,oChan2)); % PCA
                idx=find(eigenvalues/max(eigenvalues) > PCA_THRESH); % discard weak dims
                topcs=topcs(:,idx);
                b=c0(iChan,oChan2)*topcs/(topcs'*c0(oChan2,oChan2)*topcs); % matrix to project on other channels
                y(iBad,iChan)=(x(iBad,oChan2)*topcs)*b'; % projection 
            end
        end
                
                
        
        dd=y(:,iChan)-x(:,iChan); % difference from projection
        %plot([y(:,iChan),x(:,iChan)]); pause
        d=mahal(dd,dd)/thresh; % excentricity of each sample
       
        w(:,iChan)=(d<1);

    end    
    %figure(10); clf; imagescc(w);  pause
    disp(mean(w(:)))
end

x=y;

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



