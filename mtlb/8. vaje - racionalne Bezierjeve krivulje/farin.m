function q = farin(B,w)
% Opis:
%   farin vrne Farinove točke racionalne Bezierjeve krivulje, izračunane 
%   iz kontrolnih točk in uteži
%
% Definicija:
%   q = farin(B,w)
%
% Vhod:
%   B   matrika velikosti n+1 x d, v kateri vsaka vrstica predstavlja
%       eno kontrolno točko racionalne Bezierjeve krivulje stopnje n
%       v prostoru dimenzije d,
%   w   seznam uteži racionalne Bezierjeve krivulje,

%
% Izhod:
%   q   matrika velikosti n x d, v kateri i-ta vrstica predstavlja
%       i-to Farinovo točko na racionalne Bezierjeve krivulje

[n,d] = size(B); n = n-1;

q = zeros(n,d);
for i=1:n
    q(i,:) = w(i)/(w(i)+w(i+1))*B(i,:) + w(i+1)/(w(i)+w(i+1))*B(i+1,:);
end

end

