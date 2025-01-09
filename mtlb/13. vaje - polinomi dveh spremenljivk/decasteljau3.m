function d = decasteljau3(B,p)
% Opis:
%   decasteljau3 izračuna vrednost polinoma, podanega s
%   koeficienti B v točki p
%
% Definicija:
%   b = decasteljau(B,p)
%
% Vhodna podatka:
%   B   matrika velikosti n+1 x n+1, ki predstavlja koeficiente
%       polinoma dveh spremenljivk stopnje n v Bezierjevi obliki
%       (element matrike na mestu (i,j), j <= n+2-i, določa
%       koeficient polinoma z indeksom (n+2-i-j, j-1, i-1)),
%   p   baricentrične koordinate točke, v kateri računamo vrednost polinoma
%
% Izhodni podatek:
%   d   vrednost polinoma, določenega z B in p

[n, ~] = size(B); n = n-1;

P = repmat(p,n,1);
d = blossom3(B,P);

end
