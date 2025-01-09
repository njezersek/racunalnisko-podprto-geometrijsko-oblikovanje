function u = alphaparam(P,a)
% Opis:
%   alphaparam sestavi alfa parametrizacijo oziroma delitev domene na
%   podlagi točk
%
% Definicija:
%   u = alphaparam(P,a)
%
% Vhodna podatka:
%   P   matrika z m+1 vrsticami, v kateri vsaka vrstica predstavlja
%       eno točko,
%   a   parameter, ki določa alfa parametrizacijo
%
% Izhodni podatek:
%   u   seznam parametrov delitve, ki so določeni rekurzivno tako, da
%   se trenutnemu parametru iz seznama u prišteje z a potencirana
%   norma razlike zaporednih točk iz seznama P

m = length(P)-1;
u = zeros(1,m+1);
for i = 1:m
    u(i+1) = u(i) + (norm(P(i+1,:) - P(i,:)))^a;
end

end