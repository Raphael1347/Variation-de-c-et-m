function [ mk,Zinv] = DecDebyeEtZinv( Z,t,w,Zo )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%t=logspace(floor(log10(min(1./w))-1), ceil(log10(max(1./w))+1), 1000).'; 
% 
Zn=(Zo-Z)./Zo;

Zn=[real(Zn) imag(Zn)].';


A1=(t*w).^2./(1+(t*w).^2);
A2=(t*w)./(1+(t*w).^2);
A=[A1 A2].';

mk=lsqnonneg(A,Zn);

% Reconstruction du data synthetique a partir des resultats (mk) de
% l'inversion Debye
Zinv=Zo*(1-mk.'*(1-1./(1+1i*t*w)));



end

