
%{
% Simple situation with states A, B, rank(A)=9, rank(B)=9; rank(AUB)=10.

10 channels, 2 data segments.
In each each segment 9 noise sources are active, projected into
the data via distinct 9*10 mixing matrices with random coefficients.
%}
NSAMPLES=100000; % size of segment
DSR=100; % determines granularity (& minimum cluster size)
FLAGS=[]; % 'norm' or 'norm2'

noise=randn(NSAMPLES,9);
noise=nt_normcol(nt_pca(noise)); % ensure perfect decorrelation (not required)
x1=noise*randn(9,10);
x2=noise*randn(9,10);
x1=nt_normcol(x1); % normalize to remove power step
x2=nt_normcol(x2); 
x=[x1; x2]; % 

nt_cluster_jd(x,DSR);
[IDX,TODSS,SCORE]=nt_cluster_jd(x,DSR,FLAGS);
disp(['scores: ',num2str(SCORE)]);

pause;

%{
% Multiple states of rank 9, rank of concatenated=10.
%}
NSTATES=3';
x=[];
noise=randn(NSAMPLES,9);
noise=nt_normcol(nt_pca(noise)); % ensure perfect decorrelation (not required)
for iState=1:NSTATES
    x1=noise*randn(9,10);
    x1=nt_normcol(x1); % normalize to remove power step
    x=[x;x1]; % 
end

nt_cluster_jd(x,DSR);
[IDX,TODSS,SCORE]=nt_cluster_jd(x,DSR,FLAGS);
disp(['scores: ',num2str(SCORE)]);
