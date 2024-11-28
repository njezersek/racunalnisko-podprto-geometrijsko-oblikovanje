function Be = bezierelv(B,k)
% Opis:
% bezierelv izvede višanje stopnje dane Bezierjeve krivulje
%
% Definicija:
% Be = bezierelv(B,k)
%
% Vhodna podatka:
% B matrika velikosti (n+1) x d, v kateri vsaka vrstica
% predstavlja d-dimenzionalno kontrolno točko Bezierjeve
% krivulje stopnje n,
% k število, ki določa, za koliko želimo zvišati stopnjo
% dane Bezierjeve krivulje
%
% Izhodni podatek:
% Be matrika velikosti (n+k+1) x d, v kateri vsaka vrstica
% predstavlja d-dimenzionalno kontrolno točko Bezierjeve
% krvulje stopnje n+k, ki ustreza dani Bezierjevi krivulji

Be = B;
for j = 1:k
    [n,m] = size(Be);
    Bn = zeros(n+1, m);
    
    Bn(1,:) = Be(1,:);
    for i = 1:n-1
        Bn(i+1,:) = (1 - i/n) .* Be(i+1,:) + (i / n) .* Be(i, :);
    end
    Bn(n+1, :) = Be(n, :);
    
    Be = Bn;
end
