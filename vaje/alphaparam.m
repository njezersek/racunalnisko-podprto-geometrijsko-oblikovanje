function u = alphaparam(P,a)
% Opis:
% alphaparam sestavi alfa parametrizacijo oziroma delitev domene na
% podlagi točk
%
% Definicija:
% u = alphaparam(P,a)
%
% Vhodna podatka:
% P matrika z m+1 vrsticami, v kateri vsaka vrstica predstavlja
% eno točko,
% a parameter, ki določa alfa parametrizacijo
%
% Izhodni podatek:
% u seznam parametrov delitve, ki so določeni rekurzivno tako, da
% se trenutnemu parametru iz seznama u prišteje z a potencirana
% norma razlike zaporednih točk iz seznama P

[n, m] = size(P);

u = zeros(1, n);

u(1) = 0;
for i = 2:n
    u(i) = u(i-1) + norm(P(i,:) - P(i-1,:))^a;
end

end