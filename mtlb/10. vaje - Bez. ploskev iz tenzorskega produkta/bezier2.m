function [bx,by,bz] = bezier2(Bx,By,Bz,u,v)
% Opis:
%   bezier2 vrne točke na Bezierjevi ploskvi iz tenzorskega produkta
%
% Definicija:
%   [bx,by,bz] = bezier2(Bx,By,Bz,u,v)
%
% Vhodni podatki:
%   Bx, By, Bz  matrike velikosti n+1 x m+1, ki predstavljajo
%               koordinate kontrolnih točk,
%   u, v        vrstici dolžine M in N, ki predstavljata parametre
%               v smereh u in v
%
% Izhodni podatki:
%   bx, by, bz  matrike velikosti N x M, ki predstavljajo točke na
%               Bezierjevi ploskvi: [bx(J,I) by(J,I) bz(J,I)] je
%               točka pri parametrih u(I) in v(J)

[n,m] = size(Bx); n=n-1; m=m-1;
M = length(u); N = length(v);

[bx,by,bz] = deal(zeros(N,M));

for i=1:N
    kontrolne_x = [];
    kontrolne_y = [];
    kontrolne_z = [];
    for j=1:m+1
        kontrolne_x(j) = bezier(Bx(:,j),v(i));
        kontrolne_y(j) = bezier(By(:,j),v(i));
        kontrolne_z(j) = bezier(Bz(:,j),v(i));
    end
    bx(i,:) = bezier(kontrolne_x',u);
    by(i,:) = bezier(kontrolne_y',u);
    bz(i,:) = bezier(kontrolne_z',u);
end

end