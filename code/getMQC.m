function mq=getMQC(rho);
%function [mq nqc]=getMQC(rho);
% this gives the MQC intensities of a given density matrix.
% (nqc corresponds to the coherence order)
%
% Paola 2005
N=log2(length(rho));
MM=mqc(N);
mq=zeros(1,1+max(max(MM)));
%nqc=zeros(1,1+max(max(MM)));
rho=rho(:);

for k=0:max(max(MM))
  mq(k+1)=sum(abs(rho(find(MM==k))).^2);
end

mq=mq/2^(2*N);

