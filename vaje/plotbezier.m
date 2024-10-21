function p = plotbezier(B,t,c)
% Opis:
% plotbezier nariše Bezierjevo krivuljo za dane kontrolne točke in
% seznam parametrov
%
% Definicija:
% p = plotbezier(B,t,c)
%
% Vhodni podatki:
% B matrika velikosti n+1 x d, ki predstavlja kontrolne točke
% Bezierjeve krivulje stopnje n v d-dimenzionalnem prostoru,
% t seznam parametrov dolžine k, pri katerih računamo vrednost
% Bezierjeve krivulje,
% c opcijski parameter, ki določa barvo krivulje
%
% Izhodni podatek:
% p grafični objekt, ki določa krivuljo

d = size(B,2);
b = bezier(B, t);

if d == 1
    n = size(B,1)-1;
    x = (0:n)'/n;
    hold on
    p_ctl = plot(x, B);
    p = plot(t, b);
end
if d == 2
    hold on
    p_ctl = plot(B(:,1), B(:,2));
    p = plot(b(:,1), b(:,2));
end
if d == 3
    hold on
    p_ctl = plot3(B(:,1), B(:,2), B(:,3));
    p = plot3(b(:,1), b(:,2), b(:,3));
end

set(p_ctl, "Color", "black");
set(p_ctl, "Marker", "x");
set(p, "Color", c);
set(p, "LineWidth", 1.5);

end