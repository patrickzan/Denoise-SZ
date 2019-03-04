function [C,A,score]=nt_cluster1D(x);
%[C,A,energy]=nt_cluster1D(x) - cluster 1D data into 2 clusters
%
%  x: column vector or matrix of data to cluster
%
%  C: centroid pairs (one pair per column)
%  A: ownership matrix
%  score: energy/total energy, for each column

maxiter=100; 

if nargin<1; error('!'); end
if size(x,1)<2; error('too small to cluster'); end

x0=x; % save unsorted
x=sort(x);

% initialize with min,max
C1=min(x);
C2=max(x);
last_energy=inf;
energy=0;
total_energy=sum(nt_demean(x).^2);

for k=1:maxiter
    
    % cumulative energy of samples minus closest centroid 
    energy_vector= cumsum( bsxfun(@minus, x, C1).^2 )  - cumsum( bsxfun(@minus,x,C2).^2 ) ;
%     B=10;
%     energy_vector= bsxfun(@rdivide, cumsum( bsxfun(@minus, x, C1).^2 ), B+(1:size(x,1))') - bsxfun(@rdivide, cumsum( bsxfun(@minus,x,C2).^2 ), (B+(1:size(x,1))')) ;
    
    % index of lowest energy in each column
    [~,idx]=min(energy_vector,[],1);    
    
    % new centroids
    C1=mean(x(1:idx,:),1); C2=mean(x(idx+1:end,:),1);

    
    % new energy
    energy = sum(bsxfun(@minus, x(1:idx,:),C1).^2) + sum(bsxfun(@minus,x(idx+1:end,:),C2).^2); 
    score=energy./total_energy;
    if last_energy <= energy; break; end
    last_energy=energy;

    %disp(num2str([k, last_energy,energy, score, idx]))    
end

% cluster ownership labels
A=ones(size(x));

for k=1:size(x,2); y(k)=x(idx(k),k); end

thresh=bsxfun(@times,ones(size(x)),y);
A(find(x0>thresh))=2;

C=[C1;C2];

if nargout==0; clear C ; end

