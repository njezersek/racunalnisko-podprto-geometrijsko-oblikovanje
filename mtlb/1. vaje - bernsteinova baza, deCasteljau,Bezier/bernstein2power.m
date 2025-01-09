function p = bernstein2power(b)
% Opis:
%   bernstein2power pretvori polinom, predstavljen s koeficienti v
%   Bernsteinovi bazi, v polinom, predstavljen v potenčni bazi
%
% Definicija:
%   b = bernstein2power(b)
%
% Vhodni podatek:
%   p   seznam koeficientov dolžine n+1, ki po vrsti pripadajo razvoju
%       polinoma stopnje n v Bernsteinovi bazi od B_0^n do B_n^n
%
% Izhodni podatek:
%   b   seznam koeficientov dolžine n+1, ki po vrsti padajoče pripadajo 
%       razvoju polinoma stopnje n v potenčni bazi

n = length(b) - 1;
p = zeros(1,n+1);

for i = 0:n
    for j = i:n
        p(n-j+1) = p(n-j+1) + (-1)^(i+j) * nchoosek(n,j)*nchoosek(j,i) * b(i+1);
    end
end

end
