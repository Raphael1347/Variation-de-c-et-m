function Z = ColeCole( Zo,m,c,tau,w,n )
% UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


a1=(1i*tau*w);

a2=1-1./(1+a1.^c);     %[c w]
a3=repmat(a2,[1,1,n]); %[c w -]

a4=m.*a3;              %[c w m]

P1=a4(1:n,:,:);
P2=a4(n+1:end,:,:);


Z=Zo*(1-(P1+P2));   %[c w m]





