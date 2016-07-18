clear
clc
close all

% Frequences des mesures synthetiques
f=logspace(-3,4,50);
w=f.*(2*pi); % frequences angulaires associees
t=logspace(floor(log10(min(1./w))-1), ceil(log10(max(1./w))+1), 1000).';
% Parametres pour la modelisation ColeCole
% Modifier la fonction ColeCole pour prendre les arguments sous cette forme
Zo = 1000;
tau = [10.^-2, 10.^-5];

n=10; %nombre de parmètre c et m à essayer
mv=linspace(0.1,1,n);
cv=linspace(0.1,1,n);

for i=1:n;
    
    m = [mv(i),1];
    
    for j=1:n;
        
        c = [cv(j), 1];
        
        
        % Aller chercher Z (complexe) de la fonction ColeCole
        % Il serait bien que la fonction ColeCole prenne comme arguments:
        % w, Z0, m, c, tau
        Z = ColeCole(Zo,m,c,tau,w);
        
        
        
        
        
        [ mk,Zinv] = DecDebyeEtZinv( Z,t,w,Zo );
        
        mtot(i,j)=sum(mk);
        MeanTau(i,j)=exp(sum(mk.*log(t))./sum(mk));
        
        for pos=1:numel(mk)-1
            mku(pos+1)=mk(pos)+mk(pos+1);
        end
        
        mku=mku./mku(end);
        Tau60=t(find(mku>=0.6,1,'first'));
        Tau10=t(find(mku>=0.1,1,'first'));
        
        Ut(i,j)=Tau60/Tau10;
        
        RMSE(i,j) = sqrt(mean(abs(Z-Zinv)).^2);
        
%         fig=figure('name',sprintf('m=%.*f et 1, c=%.*f et 1',mv(i),cv(j)),'numbertitle','off');
%         
%         bar(t,mk)
%         set(gca,'XScale','log')
%         xlabel('Tau')
%         ylabel('mk')
%         
    end
end


figure(1)

[Mv,Cv]=meshgrid(mv,cv);

surfc(Mv,Cv,RMSE)
view(2)

% imagesc(cv,mv,RMSE)


