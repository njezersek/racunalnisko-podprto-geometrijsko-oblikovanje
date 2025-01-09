function p = plotrbezier(B,w,t,c)
% Opis:
%   plotrbezier nariše racionalno Bezierjevo krivuljo za dane kontrolne 
%   točke, uteži in seznam parametrov t
%
% Definicija:
%   p = plotbezier(B,w,t,c)
%
% Vhodni podatki:
%   B   matrika velikosti n+1 x d, ki predstavlja kontrolne točke
%       Bezierjeve krivulje stopnje n v d-dimenzionalnem prostoru,
%   w   seznam uteži racionalne Bezierjeve krivulje,
%   t   seznam parametrov dolžine k, pri katerih računamo vrednost
%       Bezierjeve krivulje
%   c   barva krivulje (opcijski parameter)
%
% Izhodni podatek:
%   p   grafični objekt, ki določa krivuljo

[n,d] = size(B); n = n-1;

% izračunamo Farinove točke:
q = farin(B,w);

% če barva krivulje ni podana, jo nastavimo na modro
if nargin == 3
    c = 'b';
end

b = rbezier(B,w,t);
hold on
if d == 2
    % kontrolni poligon:
    p = plot(B(:,1), B(:,2),'o-',Color='k',MarkerEdgeColor=c,MarkerFaceColor=c);
    % krivulja:
    plot(b(:,1), b(:,2),c);
    % Farinove točke:
    scatter(q(:,1),q(:,2),MarkerEdgeColor='r');
elseif d == 3
    x = b(:,1); y = b(:,2); z = b(:,3);
    xv = linspace(min(x), max(x), 20);
    yv = linspace(min(y), max(y), 20);
    [X,Y] = meshgrid(xv, yv);
    Z = griddata(x,y,z,X,Y,'v4');
    mesh(X,Y,Z)
    axis tight; hold on
    p = plot3(x,y,z,'.','MarkerSize',15);
    scatter3(q(:,1),q(:,2),q(:,3));
end

end