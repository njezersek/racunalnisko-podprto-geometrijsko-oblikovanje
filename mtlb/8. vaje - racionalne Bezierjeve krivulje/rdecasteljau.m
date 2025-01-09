function D = rdecasteljau(b,w,t)
% Opis:
%   rdecasteljau vrne shemo de Casteljaujevega postopka za dan seznam
%   koordinat b pri danem parametru t
%
% Definicija:
%   D = rdecasteljau(b,w,t)
%
% Vhodni podatki:
%   b   n+1 x 1 seznam koordinat kontrolnih točk Bezierjeve krivulje
%       stopnje n,
%   w   n+1 x 1 seznam uteži racionalne Bezierjeve krivulje,
%   t   parameter, pri katerem računamo koordinato Bezierjeve
%       krivulje
%
% Izhodni podatek:
%   D   tabela velikosti n+1 x n+1, ki predstavlja de Casteljaujevo
%       shemo za koordinate b in uteži w pri parametru t (element 
%       na mestu (1,n+1) je koordinata Bezierjeve krivulje pri 
%       parametru t, elementi na mestih (i,j) za i > n-j+2 so NaN)

n = length(b) - 1;
D = [b NaN(n+1,n)];
W = [w NaN(n+1,n)];

for i=1:n % stolpci - 1
    for j=0:n-i % vrstice - 1
        W(j+1,i+1) = (1-t)*W(j+1,i) + t*W(j+2,i);
        D(j+1,i+1) = (1-t)*W(j+1,i)/W(j+1,i+1)*D(j+1,i) + t*W(j+2,i)/W(j+1,i+1)*D(j+2,i);
    end
end


end