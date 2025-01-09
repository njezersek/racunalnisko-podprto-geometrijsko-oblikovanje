function B = bezierarc(fi,m)
% Opis:
%   bezierarc izračuna kontrole točke Bezierjeve krivulje, ki
%   predstavlja interpolacijo krožnega loka po izbrani metodi
%
% Definicija:
%   B = bezierarc(fi,m)
%
% Vhodna podatka:
%   fi  kot, ki določa krozni lok v območju (-fi,fi),
%   m   metoda interpolacije:
%       1 = kvadratični G1 interpolant,
%       2 = kubični C1 interpolant,
%       3 = kubični G1 interpolant s C0 interpolacijo v srednji točki
%
% Izhodni podatek:
%   B   tabela velikosti 3 x 2 ali 4 x 2, v kateri vsaka vrstica
%       predstavlja kontrolno točko Bezierjeve krivulje

if m == 1
    B = zeros(3,2);
    B(1,:) = [cos(fi) -sin(fi)];
    B(2,:) = [1/cos(fi), 0];
    B(3,:) = [cos(fi) sin(fi)];
elseif m == 2
    B = zeros(4,2);
    B(1,:) = [cos(fi) -sin(fi)];
    B(2,:) = [(cos(fi)+1/3*sin(fi)), (-sin(fi)+1/3*cos(fi))];
    B(3,:) = [(cos(fi)+1/3*sin(fi)), (sin(fi)-1/3*cos(fi))];
    B(4,:) = [cos(fi) sin(fi)];
elseif m == 3
    B = zeros(4,2);
    B(1,:) = [cos(fi) -sin(fi)];
    B(2,:) = 1/3*[(4-cos(fi)), (4*cot(fi)-4/sin(fi)+sin(fi))];
    B(3,:) = 1/3*[(4-cos(fi)), (-4*cot(fi)+4/sin(fi)-sin(fi))];
    B(4,:) = [cos(fi) sin(fi)];
end

end
