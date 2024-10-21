function [db, dB] = bezierder(B,r,t)
% Opis:
% bezierder vrne točke na krivulji, ki predstavlja odvod dane
% Bezierjeve krivulje
%
% Definicija:
% db = bezierder(B,r,t)
%
% Vhodni podatki:
% B matrika kontrolnih točk Bezierjeve krivulje, v kateri vsaka
% vrstica predstavlja eno kontrolno točko,
% r stopnja odvoda, ki ga računamo,
% t seznam parameterov, pri katerih računamo odvod
%
% Izhodni podatek:
% db matrika, v kateri vsaka vrstica predstavlja točko r-tega
% odvoda pri istoležnem parametru iz seznama t

[N, d] = size(B);
n = N-1;
k = length(t);

% izračunaj točke glede na parameter t
db = zeros(k,d);
for j = 1:k
    for i = 1:d
        D = decasteljau(B(:,i),t(j));
        vrstica = deljenadif(r, D, 1);
        db(j, i) = factorial(n)/factorial(n-r)*vrstica(n-r+1);
    end
end

% izračunaj kontrolne točke
dB = zeros(n-r+1, d);
for i = 1:n-r+1
    dB(i,:) = factorial(n)/factorial(n-r)*deljenadif(r, B, i);
end

end
