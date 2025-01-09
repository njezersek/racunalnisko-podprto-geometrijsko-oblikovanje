function u = vectorbary(P,T)
% Opis:
%   vectorbary vrne baricentrične koordinate vektorja P-0 glede na dan
%   trikotnik T
%
% Definicija: 
%   u = vectorbary(P,T)
%
% Vhodna podatka:
%   P   vrstica, ki predstavlja kartezične koordinate točke v ravnini,
%   T   matrika s tremi vrsticami, v kateri vsaka vrstica predstavlja
%       kartezične koordinate oglišč trikotnika
%
% Izhodni podatek:
%   u   vrstica dolžine 3, ki predstavlja baricentrične koordinate
%       vektorja P-0, glede na trikotnik, podan s T

u = pointbary(P,T)-pointbary([0 0],T);

end

