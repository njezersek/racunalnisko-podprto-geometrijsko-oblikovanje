function p = plotbezier2(Bx,By,Bz,u,v)
% Opis:
%   plotbezier2 nariše Bezierjevo ploskev iz tenzorskega produkta
%
% Definicija:
%   p = plotbezier2(Bx,By,Bz,u,v)
%
% Vhodni podatki:
%   Bx, By, Bz  matrike velikosti n+1 x m+1, ki predstavljajo
%               koordinate kontrolnih točk,
%   u, v        vrstici, ki predstavljata parametre v smereh u in v
%
% Izhodni podatek:
%   p   grafični objekt, ki določa ploskev

[bx,by,bz] = bezier2(Bx,By,Bz,u,v);
hold on;
mesh(Bx,By,Bz,'EdgeColor','k','FaceAlpha',0);
p = surf(bx,by,bz);

end
