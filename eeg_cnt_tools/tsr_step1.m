function [cref,cxref,r]=tsr_step1(x,ref,shifts,wx,wref,keep,thresh)
% [cref,cxref,r]=tsr_step1(x,ref,shifts,wx,wref,keep,thresh)
% calculate the auto & cross correlations
%
%  cref:  autocorrelation of reference channels
%  cxref: crosscorrelation between reference channels and data channels
%  r: channel means (removed by tsr)
% 
%  x: data to denoise (time * channels * trials)
%  ref: reference (time * channels * trials)
%  shifts: array of shifts to apply to ref (default: [0])
%  wx: weights to apply to x (time * 1 * trials);
%  wref: weights to apply to ref (time * 1 * trials);
%  keep: number of shifted-ref PCs to retain (default: all)
%  thresh: ignore shifted-ref PCs smaller than thresh (default: 10.^-12)

% Copyright 2007, 2008 Alain de Cheveigne
% modified Nai Ding 04/27/09
% Nai Ding added the outlier rejection part 07/02/09

% See: 
% de Cheveign\'e, A. and Simon, J. Z. (2007). "Denoising based on
% Time-Shift PCA." Journal of Neuroscience Methods 165: 297-305.
%
% The basic idea is to project the signal X on a basis formed by the
% orthogonalized time-shifted REF, and remove the projection. Supposing REF
% gives a good observation of the noise that contaminates X, the noise is
% removed. By allowing time shifts, the algorithm finds the optimal FIR filter 
% to apply to REF so as to compensate for any convolutional mismatch
% between X and REF.
% 
% Implementation issues:
% - Data are often available as an array of epochs. This implementation
% caters for 3D data (time * channnels * trials);
% - It is important to deemphasize high amplitude artifacts and glitches
% so that they do not dominate the solution.  This implementation uses
% weighted covariance and means.
% - Processing assumes zero-means data. Means are calculated with weights.
% - The implementation tries to be efficent and minimize memory requirements
% so as to handle large data sets.
%
% Larger data sets (disk based) could be handled by performing mean and
% covariance calculations block-by-block, in several passes.


if nargin<2; error('too few arguments'); end
if nargin<3 || isempty(shifts); shifts=0; end
if nargin<4; wx=[]; end
if nargin<5; wref=[]; end
if nargin<6 || isempty(keep); keep=[]; end
if nargin<7 || isempty(thresh); thresh=10.^-20; end

% check argument values for sanity
if size(x,1)~=size(ref,1); error('X and REF should have same nrows'); end
if size(x,3)~=size(ref,3); error('X and REF should have same npages'); end
% if ~isempty(wx) && size(x,3)~=size(wx,3); error('X and WX should have same npages'); end
% if ~isempty(wx) && size(x,1)~=size(wx,1); error('X and WX should have same nrows'); end
% if ~isempty(wref) && size(ref,1)~=size(wref,1); error('REF and WREF should have same nrows'); end
% if ~isempty(wref) && size(ref,3)~=size(wref,3); error('REF and WREF should have same npages'); end
if max(shifts)-min(0,min(shifts)) >= size(x,1); error('X has too few samples to support SHIFTS'); end
% if ~isempty(wx) && size(wx,2)~=1; error('wx should have ncols=1'); end
% if ~isempty(wref) && size(wref,2)~=1; error('wref should have ncols=1'); end
% if ~isempty(wx) && sum(wx(:))==0; error('weights on x are all zero!'); end
% if ~isempty(wref) && sum(wref(:))==0; error('weights on ref are all zero!'); end

%remove outliers
for notrial=1:size(ref,3)
    %data for each trial
    xnow=x(:,:,notrial);
    refnow=ref(:,:,notrial);
    timecourse=sum(abs(xnow),2)+sum(abs(refnow),2);
    %badtime=find(timecourse>1e5);
    badtime=find(timecourse>quantile(timecourse,0.98));
    x(badtime,:,notrial)=0;
    ref(badtime,:,notrial)=0;
end

% We need to adjust x and ref to ensure that shifts are non-negative.  If
% some values of shifts are negative, we truncate x and increment shifts.

% adjust to make shifts non-negative
n0=size(x,1);
offset1=max(0,-min(shifts));
idx=1+offset1:size(x,1);
x=x(idx,:,:);
if ~isempty(wx); wx=wx(idx,:,:); end
ref=ref(1:end-offset1,:,:);
if ~isempty(wref); wref=wref(1:end-offset1,:,:); end
shifts=shifts+offset1;                    % shifts are now positive

% adjust size of x
offset2=max(0,max(shifts)); 
idx=1: size(x,1)-offset2; 
x=x(idx,:,:);                           % part of x that overlaps with time-shifted refs
% %ref=ref(1:end-offset1,:,:);
if ~isempty(wx); wx=wx(idx,:,:); end

[mx,nx,ox]=size(x);
[mref,nref,oref]=size(ref);

% The following part is removed
% since w=1 costs more time than w=[];
% %Place Two  Nai 04/21/09  Old
% % consolidate weights into single weight matrix
% w=zeros([mx,1,oref]);
% if isempty(wx) && isempty(wref)
%     w(1:mx,:,:)=1;
% elseif isempty(wref);
%     w(:,:,:)=wx(:,:,:);
% elseif isempty(wx)
%     for k=1:ox
%         wr=wref(:,:,k);
%         wr=multishift(wr,shifts);
%         wr=min(wr,[],2);
%         w(:,:,k)=wr;
%     end;
% else;
%     for k=1:ox
%         wr=wref(:,:,k);
%         wr=multishift(wr,shifts);
%         wr=min(wr,[],2);
%         wr=min(wr,wx(1:size(wr,1),:,k));
%         w(:,:,k)=wr;
%     end
% end
% wx=w;
% wref=zeros(mref,1,oref);
% wref(idx,:,:)=w;
% %Place Two  Nai 04/21/09  Old

% remove weighted means
[x,mn1]=demean(x,wx);
ref=demean(ref,wref);

% The following part is replaced with matrix inversion
% %Place One  Nai 04/21/09  Old
% % equalize power of ref chans, then equalize power of ref PCs
% ref=normcol(ref,wref);
% ref=tspca(ref,0,[],10^-6);
% ref=normcol(ref,wref);
% 
% % covariances and cross covariance with time-shifted refs
% [cref,twcref]=tscov(ref,shifts,wref);
% [cxref,twcxref]=tsxcov(x,ref,shifts,wx);
% 
% % regression matrix of x on time-shifted refs
% r=regcov(cxref/twcxref,cref/twcref,keep,thresh);
% 
% %r=r*0.765;
% %Place One  Nai 04/21/09  Old
%Place One  Nai 04/21/09    New
%find bad samples
% badtime1=find(mean(abs(ref),2)>5e2);
% badtime2=find(mean(abs(x),2)>1e3);
% badtime=union(badtime1,badtime2);
% bad_pos_1=mod(badtime-1,size(x,1))+1;
% bad_pos_2=floor(badtime/size(x,1));
% plot(ref(bad_pos_1,:,bad_pos_2))


[cref]=tscov(ref,shifts);
[cxref]=tsxcov(x,ref,shifts); 
% r=inv(cref)*cxref';
%Place One  Nai 04/21/09    New
