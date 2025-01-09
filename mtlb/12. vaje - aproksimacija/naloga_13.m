addpath('..\1. vaje - bernsteinova baza, deCasteljau,Bezier','..\10. vaje - Bez. ploskev iz tenzorskega produkta')

% 13. naloga:
% Aproksimirali bomo funkcijo 'peaks'.
% Tabeli X in Y sta velikosti 50 × 50 in določata kartezične koordinate 
% enakomerno razporejenih točk v kvadratu [−3, 3] × [−3, 3]. Tabela Z je 
% prav tako velikosti 50 × 50, vsak njen element pa ustreza vrednosti 
% funkcije v točki, ki je določena z istoležnima elementoma iz tabel X in 
% Y. Vrednostim iz tabele Z priredimo parametre (uk, vk), k = 1,2,...,2500, 
% ki jih dobimo z linearno transformacijo kartezičnih koordinat s kvadrata 
% [−3, 3] × [−3, 3] na kvadrat [0, 1] × [0, 1]. Na ta način je določena 
% tabela podatkov P velikosti 2500 × 3.

[X,Y,Z] = peaks(50);
P = [(X(:)+3)/6 (Y(:)+3)/6 Z(:)];

% izračun kontrolnih točk Bezierjeve ploskve
m = 6; n = 5;
[Bx,By] = meshgrid(linspace(-3,3,m+1),linspace(-3,3,n+1));
Bz = lsqbezier2(m,n,P);

% izračun vrednosti v točkah mreže
[u,v] = deal(linspace(0,1,7));
[~,~,bz] = bezier2(Bx,By,Bz,u,v);

% maksimalna absolutna napaka v točkah mreže
[u,v] = deal(linspace(0,1,50));
[~,~,bz] = bezier2(Bx,By,Bz,u,v);
disp('maksimalna absolutna napaka v točkah mreže za m = 6, n = 5 :')
norm(bz(:)-P(:,3),Inf)

% izračun za višjo stopnjo
m = 9; n = 10;
[Bx,By] = meshgrid(linspace(-3,3,m+1),linspace(-3,3,n+1));
Bz = lsqbezier2(m,n,P);
[bx,by,bz] = bezier2(Bx,By,Bz,u,v);
disp('maksimalna absolutna napaka v točkah mreže za m = 9, n = 10 :')
norm(bz(:)-P(:,3),Inf)

figure
surf(X,Y,Z);
figure
surf(bx,by,bz);