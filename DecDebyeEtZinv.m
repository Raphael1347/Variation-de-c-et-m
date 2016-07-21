function [ mk,Zinv] = DecDebyeEtZinv( Z,t,w,Zo,n,numfw)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


Zn=(Zo-Z)./Zo;

Zn=[real(Zn) imag(Zn)];  %[c w m]
Zn=permute(Zn,[2,1,3]);  %[w c m]


    
A1=(t*w).^2./(1+(t*w).^2);
A2=(t*w)./(1+(t*w).^2);
A=[A1 A2].';             %[w T]

clear A1 A2


for prof=1:n;
for col=1:n;

mk(:,col,prof)=lsqnonneg(A,Zn(:,col,prof)); %[T c m]
end
end


% Reconstruction du data synthetique a partir des resultats (mk) de
% l'inversion Debye

Z12=repmat((1-1./(1+1i*t*w)),[1,1,n,n]);%[T w - -]

Z12=permute(Z12,[2 3 4 1]); %[w - - T]

mki=repmat(mk,[1,1,1,numfw]); %[T c m -]
mki=permute(mki,[4 2 3 1]);    %[- c m T]

Zinv=Zo*(1-sum(mki.*Z12,4)); %[w c m]
