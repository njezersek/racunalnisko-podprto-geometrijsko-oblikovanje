function p = plotbezier(B,t,c)
% Opis:
%   plotbezier nariše Bezierjevo krivuljo za dane kontrolne točke in
%   seznam parametrov
%
% Definicija:
%   p = plotbezier(B,t,c)
%
% Vhodni podatki:
%   B   matrika velikosti n+1 x d, ki predstavlja kontrolne točke
%       Bezierjeve krivulje stopnje n v d-dimenzionalnem prostoru,
%   t   seznam parametrov dolžine k, pri katerih računamo vrednost
%       Bezierjeve krivulje.
%   c   opcijski parameter, ki določa barvo krivulje
%
% Izhodni podatek:
%   p   grafični objekt, ki določa krivuljo

if nargin == 2
    c = 'b';
end

b = bezier(B,t);
d = length(B(1,:));
hold on
if d == 2
    p = plot(B(:,1), B(:,2),'o-',Color='k',MarkerEdgeColor=c,MarkerFaceColor=c);
    plot(b(:,1), b(:,2),c);
elseif d == 3
    x = b(:,1); y = b(:,2); z = b(:,3);
    xv = linspace(min(x), max(x), 20);
    yv = linspace(min(y), max(y), 20);
    [X,Y] = meshgrid(xv, yv);
    Z = griddata(x,y,z,X,Y,'v4');
    mesh(X,Y,Z)
    axis tight; hold on
    p = plot3(x,y,z,'.','MarkerSize',15);
end

end