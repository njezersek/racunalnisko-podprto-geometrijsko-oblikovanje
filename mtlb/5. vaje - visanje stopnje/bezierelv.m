function Be = bezierelv(B,k)
% Opis:
%   bezierelv izvede višanje stopnje dane Bezierjeve krivulje
%
% Definicija:
%   Be = bezierelv(B,k)
%
% Vhodna podatka:
%   B   matrika velikosti (n+1) x d, v kateri vsaka vrstica
%       predstavlja d-dimenzionalno kontrolno točko Bezierjeve
%       krivulje stopnje n,
%   k   število, ki določa, za koliko želimo zvišati stopnjo
%       dane Bezierjeve krivulje
%
% Izhodni podatek:
%   Be  matrika velikosti (n+k+1) x d, v kateri vsaka vrstica
%       predstavlja d-dimenzionalno kontrolno točko Bezierjeve
%       krvulje stopnje n+k, ki ustreza dani Bezierjevi krivulji

for i=1:k
    Be = bezierelv_1(B);
    B = Be;
end


end