function S = argyrisspline(tri,f,Df,Hf)
% Opis:
%   argyrisspline izračuna Bezierjeve predstavitve polinomov, ki
%   določajo Argyrisov zlepek nad triangulacijo
%
% Definicija:
%   S = argyrisspline(tri,f,Df,Hf)
%
% Vhodni podatki:
%   tri     objekt razreda triangulation, ki določa triangulacijo
%           domene, nad katero aproksimiramo funkcijo f,
%   f       funkcija, ki jo interpoliramo,
%   Df      gradient funkcije, ki jo interpoliramo,
%   Hf      Hessejeva matrika funkcije, ki jo interpoliramo
%
% Izhodni podatek:
%   S       celica z dolžino, ki ustreza številu trikotnikov v
%           triangulaciji, v kateri vsak element vsebuje matriko
%           velikosti 6 x 6, ki predstavlja koeficiente Argyrisovega
%           polinoma v Bezierjevi obliki

n = length(tri.ConnectivityList);
S = cell(n);
for i=1:n
    T = tri.ConnectivityList(i,:);
    T = [tri.Points(T(1),:); tri.Points(T(2),:); tri.Points(T(3),:)];
    S{i} = argyris(T,f,Df,Hf);
end

end