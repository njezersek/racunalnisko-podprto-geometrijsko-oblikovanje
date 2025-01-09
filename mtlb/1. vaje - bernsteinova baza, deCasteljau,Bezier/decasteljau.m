function D = decasteljau(b,t)
% Opis:
%   decasteljau vrne shemo de Casteljaujevega postopka za dan seznam
%   koordinat b pri danem parametru t
%
% Definicija:
%   D = decasteljau(b,t)
%
% Vhodna podatka:
%   b   seznam koordinat kontrolnih točk Bezierjeve krivulje
%       stopnje n,
%   t   parameter, pri katerem računamo koordinato Bezierjeve
%       krivulje
%
% Izhodni podatek:
%   D   tabela velikosti n+1 x n+1, ki predstavlja de Casteljaujevo
%       shemo za koordinate b pri parametru t (element na mestu
%       (1,n+1) je koordinata Bezierjeve krivulje pri parametru t,
%       elementi na mestih (i,j) za i > n-j+2 so NaN)

n = length(b) - 1;
D = [b NaN(n+1,n)];

for i=1:n % stolpci - 1
    for j=0:n-i % vrstice - 1
        D(j+1,i+1) = (1-t).*D(j+1, i) + t .* D(j + 2, i);
    end
end


end