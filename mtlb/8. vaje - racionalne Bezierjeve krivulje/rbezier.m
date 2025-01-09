function b = rbezier(B,w,t)
% Opis:
%   rbezier vrne točke na racionalni Bezierovi krivulji, izračunane z
%   de Casteljaujevim postopkom za racionalne krivulje
%
% Definicija:
%   b = rbezier(B,w,t)
%
% Vhod:
%   B   matrika velikosti n+1 x d, v kateri vsaka vrstica predstavlja
%       eno kontrolno točko racionalne Bezierjeve krivulje stopnje n
%       v prostoru dimenzije d,
%   w   seznam uteži racionalne Bezierjeve krivulje,
%   t   seznam parametrov dolžine N, za katere se računajo točke na
%       racionalni Bezierjevi krivulji
%
% Izhod:
%   b   matrika velikosti N x d, v kateri i-ta vrstica predstavlja
%       točko na racionalni Bezierjevi krivulji pri i-tem parametru
%       iz seznama t


[n,d] = size(B); n=n-1;
N = length(t);

b = zeros(N,d);
D = cell(1,d);

for i=1:N
    for j = 1:d
        D1 = rdecasteljau(B(:,j),w,t(i));
        D{j} = D1;
        b(i,j) = D1(1,end);
    end
end

end
