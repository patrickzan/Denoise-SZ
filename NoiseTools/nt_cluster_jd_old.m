function [A,todss]=nt_cluster_jd(x,dsr,flags)
%[A,todss]=nt_cluster_jd(x,dsr,flags) - cluster with joint diagonalization
%
%  A: map of cluster ownership
%  todss: result of DSS
%
%  x: data (time*channels)
%  dsr: downsample ratio for cross product series
%  flags: 'norm': give each dsr-sized slice the same weight
%
% See nt_bias_cluster, nt_cluster1D


if nargin<3 ||isempty(flags); flags=[]; end
if nargin<2; error('!'); end

if ~exist('vl_kmeans');
    disp('vl_kmeans() not found, download from http://www.vlfeat.org');
end
if ndims(x)>2; 
    error('x should be time*channels');
end

disp('entering nt_cluster_jd...');

% initial clustering
[C0,C1,A,todss]=nt_bias_cluster(x,dsr,flags);
todss2=todss(:,[1 end]); % keep only first and last components
z=nt_mmat(x,todss2);

figure(2);  clf; set(gcf, 'name','nt_cluster_jd');
subplot 511; plot(x); title('data');
subplot 512; plot(A,'.-'); title('initial cluster map');
subplot 513; plot(z(:,1)); title('DSS1');
subplot 514; plot(z(:,2)); title('DSS2');
drawnow; %pause(1);

% iterate until stable
old_A=A;
for k=1:10
    [~,~,A]=nt_bias_cluster(nt_normcol(z),dsr,flags); %cluster first & last components
    disp(['cluster sizes: ', num2str([sum(A==1), sum(A==2)])]);
    c1=nt_cov(x(find(A==1),:));
    c0=c1+nt_cov(x(find(A==2),:));
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    todss2=todss(:,[1 end]); % keep only first and last
    z=nt_mmat(x,todss2);
    %figure(3); clf; plot(double(A)); pause
    subplot 515; plot(A,'.-'); title('final cluster map');
    if old_A==A; break; end
    old_A=A;
end 

if nargout==0;
    % no output, just plot
    disp(['cluster1: ',num2str(100*numel(find(A==1))/numel(A)), '%']);

    z1=nt_mmat(x,todss(:,1));
    z2=nt_mmat(x,todss(:,end));
    z=nt_mmat(x,todss); 
    z=nt_normcol(z);
    e1=mean(z(find(A==1),:).^2);
    e2=mean(z(find(A==2),:).^2);

    figure(100); clf
    plot(x); hold on
    xx=x;
    xx(find(A==2),:)=nan;
    plot(xx,'k');
    axis tight
    title('black: cluster2');
    
    figure(101); clf
    subplot 121;
    plot(pwr1./pwr0,'.-'); xlabel('component'); ylabel('score'); title('DSS cluster A vs B');
    subplot 122;
    wsize=min(1024,size(z1,1));
    nt_spect_plot(z1,wsize,[],[],1);
    hold on
    wsize=min(1024,size(z2,1));
    nt_spect_plot(z2,wsize,[],[],1);
    xlim([0 .5])
    nt_colorlines([],[1 3]);
    legend('first','last'); legend boxoff
    hold off

    
    figure(102); clf
    subplot 211;
    plot(z1); axis tight
    title('first DSS component')
    subplot 212;
    plot(z2); axis tight
    title('last DSS component');
    
    figure(103); clf
    plot([e1',e2'], '.-'); legend('cluster A', 'cluster B'); title ('power per component');
    
    figure(104);clf
    subplot 121; imagescc(nt_cov(z(find(A==1),:))); title('cluster A'); 
    subplot 122; imagescc(nt_cov(z(find(A==2),:))); title('cluster B');
    
    if 0 
        figure(105); clf
        subplot 211;
        nt_sgram(z1,1024,32,[],1);
        title('first');
        subplot 212;
        nt_sgram(z2,1024,32,[],1);
        title('last');
    end
    clear c0 c1 A todss pwr0 pwr1
    
end

function y=norm2(x,n,ind)
[I,J]=ind2sub([n,n],ind);
for k=1:size(x,1)
    a=x(k,1:n);
    b=sqrt(a(I).*a(J));
    y(k,:)=x(k,:)./b;
end

    
    
    