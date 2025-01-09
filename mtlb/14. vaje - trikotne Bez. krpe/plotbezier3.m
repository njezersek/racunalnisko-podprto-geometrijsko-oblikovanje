function p = plotbezier3(Bx,By,Bz,d)
% Opis:
%   plotbezier nariše Bezierjevo ploskev iz tenzorskega produkta
%
% Definicija:
%   p = plotbezier3(Bx,By,Bz,d)
%
% Vhodni podatki:
%   Bx, By, Bz  matrike velikosti n+1 x n+1, ki predstavljajo
%               koordinate kontrolnih točk Bezierjeve krpe (element
%               posamezne matrike na mestu (i,j), j <= n+2-i, določa
%               koordinato kontrolne točke z indeksom
%               (n+2-i-j, j-1, i-1)),
%   d           vhodni podatek za metodo trimeshgrid
%
% Izhodni podatek:
%   p           grafični objekt, ki določa ploskev

hold on;

% Narišemo ploskev:
[TRI,U] = trimeshgrid(d);
b = bezier3(Bx,By,Bz,U);
p = trisurf(TRI, b(:,1), b(:,2), b(:,3));

% Narišemo kontrolni poligon:
[n, ~] = size(Bx); n = n - 1; % Red ploskve
controlPoints = []; % Shramba za kontrolne točke

for i = 1:n+1
    for j = 1:n+2-i
        % Zgradimo kontrolne točke v 3D
        controlPoints = [controlPoints; Bx(i,j), By(i,j), Bz(i,j)];
    end
end

% Triangulacija kontrolnih točk
[controlTRI, controlU] = trimeshgrid(n);
controlTRI(controlTRI > size(controlPoints, 1)) = []; % Odstranimo presežke
trimesh(controlTRI, controlPoints(:,1), controlPoints(:,2), controlPoints(:,3), 'EdgeColor', 'k', 'FaceAlpha',0); % Črtkani črni robovi


end


