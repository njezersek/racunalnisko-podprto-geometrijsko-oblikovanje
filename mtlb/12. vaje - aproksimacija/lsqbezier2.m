function B = lsqbezier2(m,n,P)
% Opis:
%   lsqbezier2 vrne koeficiente tenzorskega polinoma, ki po metodi
%   najmanjših kvadratov najbolje aproksimira dane podatke
%
% Definicija:
%   B = lsqbezier2(m,n,P)
%
% Vhodni podatki:
%   m,n     parametra, ki določata stopnjo polinoma,
%   P       matrika podatkov, ki v vsaki vrstici vsebuje parametra z
%           intervala [0,1] ter njima pripadajočo vrednost
%
% Izhodni podatek:
%   B       matrika velikosti n+1 x m+1 s koeficienti polinoma, ki po
%           metodi najmanjših kvadratov najbolje aproksimira podatke

bern = @(n,i,x) nchoosek(n,i) * x^i * (1-x)^(n-i);

K = size(P,1);
A = NaN(K, (m+1)*(n+1));

for i=1:K
    for j1 = 0:m
        for j2 = 0:n
            A(i,j1*(n+1)+(j2+1)) = bern(m,j1,P(i,1)) * bern(n,j2,P(i,2));
        end
    end
end

b = P(:,3);
B = A\b;
B = reshape(B,[n+1, m+1]); % iz stolpca naredi n+1 x m+1 matriko

end