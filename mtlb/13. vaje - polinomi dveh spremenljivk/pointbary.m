function u = pointbary(P,T)
% Opis:
%   pointbary vrne baricentrične koordinate točke glede na dan trikotnik T
%
% Definicija:
%   u = pointbary(P,T)
%
% Vhodna podatka:
%   P   vrstica, ki predstavlja kartezične koordinate točke v ravnini,
%   T   matrika s tremi vrsticami, v kateri vsaka vrstica predstavlja
%       kartezične koordinate oglišč trikotnika
%
% Izhodni podatek:
%   u   vrstica dolžine 3, ki predstavlja baricentrične koordinate
%       točke, podane s P, glede na trikotnik, podan s T

A = [1 1 1; T(1,:)' T(2,:)' T(3,:)'];
b = [1;P(1);P(2)];

u = A\b;
u = u';

end
