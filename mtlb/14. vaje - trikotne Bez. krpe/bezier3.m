function b = bezier3(Bx,By,Bz,U)
% Opis:
%   bezier3 izračuna točke na trikotni Bezierjevi ploskvi
%
% Definicija:
%   b = bezier3(Bx,By,Bz,U)
%
% Vhodni podatki:
%   Bx, By, Bz  matrike velikosti n+1 x n+1, ki predstavljajo
%               koordinate kontrolnih točk Bezierjeve krpe (element
%               posamezne matrike na mestu (i,j), j <= n+2-i, določa
%               koordinato kontrolne točke z indeksom
%               (n+2-i-j, j-1, i-1)),
%   U           matrika, v kateri vrstice predstavljajo baricentrične
%               koordinate točk glede na domenski trikotnik, za katere
%               računamo točke na Bezierjevi krpi
%
% Izhodni podatek:
%   b           matrika, v kateri vsaka vrstica predstavlja točko na
%               Bezierjevi krpi pri istoležnih parametrih iz matrike U

[n,~] = size(Bx); n = n-1;
[k,~] = size(U);
d = 3;

b = zeros(k,d);

for i=1:k
    dx = decasteljau3(Bx,U(i,:));
    dy = decasteljau3(By,U(i,:));
    dz = decasteljau3(Bz,U(i,:));
    b(i,:) = [dx dy dz];
end

end