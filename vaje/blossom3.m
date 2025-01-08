function b = blossom3(B, U)
% Opis:
% blossom3 izračuna razcvet polinoma dveh spremenljivk
%
% Definicija:
% b = blossom3(B, U)
%
% Vhodna podatka:
% B matrika velikosti n+1 x n+1, ki predstavlja koeficiente
% polinoma dveh spremenljivk stopnje n v Bezierjevi obliki
% (element matrike na mestu (i,j), j <= n+2-i, določa
% koeficient polinoma z indeksom (n+2-i-j, j-1, i-1)),
% U matrika velikosti n x 3, v kateri vrstice predstavljajo
% baricentrične koordinate točk glede na domenski trikotnik,
% za katere izvajamo razcvet polinoma
%
% Izhodni podatek:
% b vrednost razcveta polinoma, določenega z matriko B, v
% točkah, določenih z matriko U

% Preverimo velikost matrike B, da ugotovimo stopnjo n
[n_plus_1, ~] = size(B);
n = n_plus_1 - 1;

% Inicializiramo delovno matriko
% Delovna matrika bo vsebovala koeficiente za razcvet
work_matrix = B;

% Izvajamo razcvet glede na vrstice matrike U
for k = 1:n
    % Baricentrične koordinate trenutne točke
    u1 = U(k, 1);
    u2 = U(k, 2);
    u3 = U(k, 3);

    % Posodobimo delovno matriko za trenutni korak
    new_matrix = zeros(n + 1 - k, n + 1 - k);
    for i = 1:(n + 1 - k)
        for j = 1:(n + 1 - k - i + 1)
            % Baricentrična kombinacija trenutnih koeficientov
            new_matrix(i, j) = u1 * work_matrix(i, j) + ...
                               u2 * work_matrix(i, j + 1) + ...
                               u3 * work_matrix(i + 1, j);
        end
    end
    work_matrix = new_matrix; % Posodobimo delovno matriko
end

% Končni rezultat je edini preostali element v delovni matriki
b = work_matrix(1, 1);
end
