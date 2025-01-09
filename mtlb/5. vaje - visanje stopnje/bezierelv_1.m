function Be = bezierelv_1(B)
% Opis:
%   bezierelv_1 izvede en korak višanja stopnje dane Bezierjeve krivulje
%
% Definicija:
%   Be = bezierelv_1(B)
%
% Vhodna podatka:
%   B   matrika velikosti (n+1) x d, v kateri vsaka vrstica
%       predstavlja d-dimenzionalno kontrolno točko Bezierjeve
%       krivulje stopnje n
%
% Izhodni podatek:
%   Be  matrika velikosti (n+2) x d, v kateri vsaka vrstica
%       predstavlja d-dimenzionalno kontrolno točko Bezierjeve
%       krvulje stopnje n+1, ki ustreza dani Bezierjevi krivulji

[n,d] = size(B); n = n-1;
Be = [B(1,:); zeros(n,d); B(end,:)];
for i = 1:n
    Be(i+1,:) = (1-i/(n+1))*B(i+1,:) + i/(n+1)*B(i,:);
end

end