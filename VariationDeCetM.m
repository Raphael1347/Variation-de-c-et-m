clear
clc
close all
tic
numfw=50; %nombre d'élément dans le vecteur de fréquence
f=logspace(-3,4,numfw);
w=f.*(2*pi); % frequences angulaires associees
t=logspace(floor(log10(min(1./w))-1), ceil(log10(max(1./w))+1), 1000).';
% Parametres pour la modelisation ColeCole
% Modifier la fonction ColeCole pour prendre les arguments sous cette forme
Zo = 1000;

n=30; %nombre de parmètre c et m à essayer

tau1 = ones(n,1)*10.^-2;
tau2 = ones(n,1)*10.^-5;
tau=[tau1 ; tau2];%pour avoir un vecteur colone contenant Tau1 suivie de Tau2

m1(1,1,:)=linspace(0.1,1,n);
m2=ones(1,1,n);

m=repmat(m1,n,numfw);
m=[m ; repmat(m2,n,numfw)];%pour avoir un vecteur dans la dimantion 3
%(profondeur) contenant m1 et m2 et répété ce
%ce vecteur pour avoir une matrice 3d de la même
%hauteur que c et de la même largeur que w
%[- - m]

c1(:,1)=linspace(0.1,1,n);
c2=ones(n,1);
c=[c1;c2]*ones(1,numfw); %pour avoir un vecteur colone contenant
%c1 suivie de c2 et répété cette colonne
%pour avoir une matrice de même largeur que w


% Aller chercher Z (complexe) de la fonction ColeCole
% Il serait bien que la fonction ColeCole prenne comme arguments:
% w, Z0, m, c, tau


Z = ColeCole(Zo,m,c,tau,w,n);

clear a1 a2 a3 P1 P2



[ mk,Zinv] = DecDebyeEtZinv( Z,t,w,Zo,n,numfw);


lnt=repmat(log(t)*ones(1,n),[1 1 n]); %créer un matrice de ln(t) de même taille que
%la matrice mk où chaque colonne est identique

MeanTau=exp(sum(mk.*lnt,1)./sum(mk,1));

[rowmk,colmk]=size(mk);
mkU=mk;


for row=1:rowmk-1;
    mkU(row+1,:,:)=mkU(row,:,:)+mkU(row+1,:,:);
end


mkUmax=repmat(mkU(end,:,:),[rowmk 1 1]);
mkU=mkU./mkUmax;


for prof=1:n
    for col=1:n
        Tau60(1,col,prof)=t(find(mkU(:,col,prof)>=0.6,1,'first'));
        Tau10(1,col,prof)=t(find(mkU(:,col,prof)>=0.1,1,'first'));
    end
end
Ut=Tau60./Tau10;


Z=permute(Z,[2 1 3]); %[w c m]

RMSE = permute(sqrt(mean(abs(Z-Zinv).^2,1)),[2 3 1]);


%         fig=figure('name',sprintf('m=%.*f et 1, c=%.*f et 1',mv(i),cv(j)),'numbertitle','off');
%
%         bar(t,mk)
%         set(gca,'XScale','log')
%         xlabel('Tau')
%         ylabel('mk')
%



figure(1)

[M1,C1]=meshgrid(m1,c1);

surfc(M1,C1,RMSE)
view(2)
xlabel('c')
ylabel('m')
zlabel('RMSE')
% imagesc(c1,m2,RMSE)
toc

