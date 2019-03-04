function [topcs,eigenvalues]=nt_pcarot(cov,N)
% [topcs,eigenvalues]=pcarot(cov,N) - PCA matrix from covariance
%
%  topcs: PCA rotation matrix
%  eigenvalues: PCA eigenvalues
%  
%  cov: covariance matrix
%  N: eigs' K parameter (if absent: use eig)
%
% NoiseTools

if nargin >1 && ~isempty(N); 
    [V, S] = eigs(cov,N) ;  
else
    [V, S] = eig(cov) ;  
end

V=real(V);
S=real(S);
[eigenvalues, idx] = sort(diag(S)', 'descend') ;
topcs = V(:,idx);

