function B = argyris(T,f,Df,Hf)
% Opis:
%   argyris izračuna Bezierjeve ordinate polinoma dveh spremenljivk
%   stopnje 5, ki interpolira vrednosti, prve in druge odvode podane
%   funkcije f v ogliščih trikotnika ter odvode funkcije f v središčih
%   stranic trikotnika v smereh, pravokotnih na posamezno stranico
%
% Definicija:
%   B = argyris(T,f,Df,Hf)
%
% Vhodni podatki:
%   T   tabela velikosti 3 x 2, v kateri vsaka vrstica predstavlja
%       oglišče trikonika, nad katerim je definiran polinom,
%   f   funkcija, ki jo interpoliramo,
%   Df  gradient funkcije, ki jo interpoliramo,
%   Hf  Hessejeva matrika funkcije, ki jo interpoliramo
%
% Izhodni podatek:
%   B   matrika velikosti 6 x 6, ki predstavlja koeficiente polinoma
%       dveh spremenljivk stopnje 5 v Bezierjevi obliki (element
%       matrike na mestu (i,j), j <= 7-i, določa koeficient polinoma
%       z indeksom (7-i-j, j-1, i-1)), ki interpolira funkcijo f

V1 = T(1,:)';
V2 = T(2,:)';
V3 = T(3,:)';

B = NaN(6);
B(1,1) = f(V1(1),V1(2));
B(1,2) = f(V1(1),V1(2)) + 1/5*dot(Df(V1(1),V1(2)),V2-V1);
B(2,1) = f(V1(1),V1(2)) + 1/5*dot(Df(V1(1),V1(2)),V3-V1);
B(1,3) = f(V1(1),V1(2)) + 2/5*dot(Df(V1(1),V1(2)),V2-V1) + 1/20*dot(V2-V1,Hf(V1(1),V1(2))*(V2-V1));
B(2,2) = f(V1(1),V1(2)) + 1/5*dot(Df(V1(1),V1(2)),V2-V1) + 1/5*dot(Df(V1(1),V1(2)),V3-V1) + 1/20*dot(V3-V1,Hf(V1(1),V1(2))*(V2-V1));
B(3,1) = f(V1(1),V1(2)) + 2/5*dot(Df(V1(1),V1(2)),V3-V1) + 1/20*dot(V3-V1,Hf(V1(1),V1(2))*(V3-V1));

B(1,6) = f(V2(1),V2(2));
B(2,5) = f(V2(1),V2(2)) + 1/5*dot(Df(V2(1),V2(2)),V3-V2);
B(1,5) = f(V2(1),V2(2)) + 1/5*dot(Df(V2(1),V2(2)),V1-V2);
B(3,4) = f(V2(1),V2(2)) + 2/5*dot(Df(V2(1),V2(2)),V3-V2) + 1/20*dot(V3-V2,Hf(V2(1),V2(2))*(V3-V2));
B(2,4) = f(V2(1),V2(2)) + 1/5*dot(Df(V2(1),V2(2)),V3-V2) + 1/5*dot(Df(V2(1),V2(2)),V1-V2) + 1/20*dot(V1-V2,Hf(V2(1),V2(2))*(V3-V2));
B(1,4) = f(V2(1),V2(2)) + 2/5*dot(Df(V2(1),V2(2)),V1-V2) + 1/20*dot(V1-V2,Hf(V2(1),V2(2))*(V1-V2));

B(6,1) = f(V3(1),V3(2));
B(5,1) = f(V3(1),V3(2)) + 1/5*dot(Df(V3(1),V3(2)),V1-V3);
B(5,2) = f(V3(1),V3(2)) + 1/5*dot(Df(V3(1),V3(2)),V2-V3);
B(4,1) = f(V3(1),V3(2)) + 2/5*dot(Df(V3(1),V3(2)),V1-V3) + 1/20*dot(V1-V3,Hf(V3(1),V3(2))*(V1-V3));
B(4,2) = f(V3(1),V3(2)) + 1/5*dot(Df(V3(1),V3(2)),V1-V3) + 1/5*dot(Df(V3(1),V3(2)),V2-V3) + 1/20*dot(V2-V3,Hf(V3(1),V3(2))*(V1-V3));
B(4,3) = f(V3(1),V3(2)) + 2/5*dot(Df(V3(1),V3(2)),V2-V3) + 1/20*dot(V2-V3,Hf(V3(1),V3(2))*(V2-V3));

a = V2-V1; n12 = 1/norm(a)*[a(2); -a(1)];
d = vectorbary(n12,T);
d1 = d(1); d2 = d(2); d3 = d(3);
V12 = 1/2*(V1 + V2);
B(2,3) = 8/15*1/d3*dot(Df(V12(1),V12(2)),n12) - 1/6*(B(2,1)+4*B(2,2)+4*B(2,4)+B(2,5)) - 1/6*d1/d3*(B(1,1)+4*B(1,2)+6*B(1,3)+4*B(1,4)+B(1,5)) - 1/6*d2/d3*(B(1,6)+4*B(1,5)+6*B(1,4)+4*B(1,3)+B(1,2));

a = V3-V2; n23 = 1/norm(a)*[a(2); -a(1)];
d = vectorbary(n23,T);
d1 = d(1); d2 = d(2); d3 = d(3);
V23 = 1/2*(V2 + V3);
B(3,3) = 8/15*1/d1*dot(Df(V23(1),V23(2)),n23) - 1/6*(B(1,5)+4*B(2,4)+4*B(4,2)+B(5,1)) - 1/6*d2/d1*(B(1,6)+4*B(2,5)+6*B(3,4)+4*B(4,3)+B(5,2)) - 1/6*d3/d1*(B(6,1)+4*B(5,2)+6*B(4,3)+4*B(3,4)+B(2,5));

a = V1-V3; n31 = 1/norm(a)*[a(2); -a(1)];
d = vectorbary(n31,T);
d1 = d(1); d2 = d(2); d3 = d(3);
V31 = 1/2*(V3 + V1);
B(3,2) = 8/15*1/d2*dot(Df(V31(1),V31(2)),n31) - 1/6*(B(1,2)+4*B(2,2)+4*B(4,2)+B(5,2)) - 1/6*d1/d2*(B(1,1)+4*B(2,1)+6*B(3,1)+4*B(4,1)+B(5,1)) - 1/6*d3/d2*(B(6,1)+4*B(5,1)+6*B(4,1)+4*B(3,1)+B(2,1));



end
