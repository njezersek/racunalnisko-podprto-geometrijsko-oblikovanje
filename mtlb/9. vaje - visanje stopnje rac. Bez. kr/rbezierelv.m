function [Be,we] = rbezierelv(B,w)
% Opis:
%   rbezierelv izvede višanje stopnje racionalne Bezierjeve krivulje
%
% Definicija:
%   [Be,we] = rbezierelv(B,w)
%
% Vhodna podatka:
%   B   matrika velikosti n+1 x d, v kateri vsaka vrstica predstavlja
%       d-dimenzionalno kontrolno točko racionalne Bezierjeve krivulje
%       stopnje n,
%   w   seznam uteži racionalne Bezierjeve krivulje
%
% Izhodni podatek:
%   Be  matrika velikosti n+2 x d, v kateri vsaka vrstica predstavlja
%       d-dimenzionalno kontrolno točko racionalne Bezierjeve krvulje
%       stopnje n+1, ki ustreza dani racionalni Bezierjevi krivulji,
%   we  seznam dolžine n+2, v katerem vsak element predstavlja utež
%       racionalne Bezierjeve krvulje stopnje n+1, ki ustreza dani
%       racionalni Bezierjevi krivulji

B1 = [w.*B w];
B1e = bezierelv_1(B1); % v mapi 5
we = B1e(:,end);
Be = B1e(:,1:(end-1))./we;

end
