function Z = ColeCole( Zo,m,c,t,w )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here




for n=1:numel(m)

z(n,:)=m(n).*(1-1./(1+(1i.*w.*t(n)).^c(n)));

end



Z=Zo.*(1-sum(z));



end