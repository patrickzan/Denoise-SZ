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

y=zeros(size(x));
d=ones(size(x,1),nchan); % distance from projection, d>thresh means artifact

NITER=2; % iterations to refine c0
for k=1:NITER
    
    % c0: covariance of non-artifact part
    dd=max(d,[],2);
    c0=nt_cov(x,[],dd<=1);
    
    d0=d;
    for iChan=1:nchan

        oChan=setdiff(1:nchan,iChan); % other channels
        
        [topcs,eigenvalues]=nt_pcarot(c0(oChan,oChan)); % PCA
        PCA_THRESH=10^-6;
        idx=find(eigenvalues/max(eigenvalues) > PCA_THRESH); % discard weak dims
        topcs=topcs(:,idx);
        b=c0(iChan,oChan)*topcs/(topcs'*c0(oChan,oChan)*topcs); % matrix to project on other channels
        y(:,iChan)=(x(:,oChan)*topcs)*b'; % projection 

        y0=y;
        
        % fix to exclude the most serious competing channels from projection
        competitors=zeros(size(d0(:,oChan)));
        iCompetitors=d0(:,oChan)> 1 ... % bad
            & d0(:,oChan) > repmat(d0(:,iChan),1,nchan-1); % and worse than this one
        competitors(find(iCompetitors))=d0(find(iCompetitors));
        idx=find(~min(competitors,[],2)); % samples with no competitors
        [~,worst]=max(competitors,[],2); % which competitor dominates?
        worst(idx)=0; % set to zero if there is none
        for k=1:numel(oChan);
            
            iReplace=find(k==worst);
            if iReplace
                ooChan=setdiff(1:nchan,[iChan,oChan(k)]); % channels other than this and worst
                [topcs,eigenvalues]=nt_pcarot(c0(ooChan,ooChan)); % PCA
                PCA_THRESH=10^-6;
                idx=find(eigenvalues/max(eigenvalues) > PCA_THRESH); % discard weak dims
                topcs=topcs(:,idx);
                b=c0(iChan,ooChan)*topcs/(topcs'*c0(ooChan,ooChan)*topcs);
                y(iReplace,iChan)=(x(iReplace,ooChan)*topcs)*b'; % projection
            end
        end
        
        
        ddd=y(:,iChan)-x(:,iChan); % difference from projection
        
        d(:,iChan)=mahal(ddd,ddd)/thresh; % excentricity of each sample and each channel 

        figure(1); clf
        FOCUS=861600+(1:500); subplot 311; plot(x(FOCUS,:)); title('raw');
        subplot 312; plot([x(FOCUS,iChan), y0(FOCUS,iChan), y(FOCUS,iChan)]); title (num2str(iChan)); legend('raw','proj','fixed');
        subplot 313; plot(d(FOCUS,:)); pause

    end    
    %disp(mean(w))
end

dd=max(d,[],2);
%x=nt_demean(x,dd<1);
c0=nt_cov(x,[],dd<1);


%{
For each channel, find the part for which it dominates the artifact.
Replace that part based on projection on the other channels. 
%}

idx=find(dd>1);
xx=x(idx,:);
xx=abs(nt_normcol(xx));

for iChan=1:nchan
    
%     oChan=setdiff(1:nchan,iChan);
%     b=c0(iChan,oChan)/c0(oChan,oChan);    % projection matrix
%     idx2=find(xx(:,iChan)>max(xx(:,oChan),[],2)); % find sections where this channel is largest
%     idx2=find(xx(:,iChan)>max(xx(:,oChan),[],2)); % find sections where this channel is largest
%     x(idx(idx2),iChan)=x(idx(idx2),oChan)*b'; % replace by projection

    oChan=setdiff(1:nchan,iChan);
    b=c0(iChan,oChan)/c0(oChan,oChan);    % projection matrix
    idx2=find(xx(:,iChan)>max(xx(:,oChan),[],2)); % find sections where this channel is largest
    
    x(idx(idx2),iChan)=x(idx(idx2),oChan)*b'; % replace by projection
    
    idx2=find(d(:,iChan)>1 & d(:,iChan)>max(d(:,oChan),[],2));
    %x(idx2,iChan)=x(idx2,oChan)*b';

end

%{
To do:
- avoid problems with rank-deficient data
- exclude borderline channels from projection
%}

x=bsxfun(@times,x,nn);
x=bsxfun(@plus,x,mn);

x=nt_fold(x,nsample);
d=nt_fold(d,nsample);



