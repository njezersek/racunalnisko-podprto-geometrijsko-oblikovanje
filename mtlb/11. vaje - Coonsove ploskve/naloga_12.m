addpath('..\1. vaje - bernsteinova baza, deCasteljau,Bezier', '..\2. vaje', '..\4. vaje - subdivizija', '..\5. vaje - visanje stopnje','..\6. vaje - parametrizacija sest. krivulje','..\8. vaje - racionalne Bezierjeve krivulje','..\10. vaje - Bez. ploskev iz tenzorskega produkta')


% 12. NALOGA:

% 1)

% kontrole točke robnih krivulj
[Bx,By,Bz] = deal(zeros(9,7));

phi = linspace(0,pi,7);
psi = linspace(0,pi,9);

Bx(1,:) = 1+sin(phi);
By(1,:) = -cos(phi);

Bx(9,:) = -1-sin(phi);
By(9,:) = -cos(phi);

Bx(:,1) = cos(psi);
By(:,1) = -1;

Bz(:,1) = sin(psi);
Bx(:,7) = cos(psi);

By(:,7) = 1;
Bz(:,7) = sin(psi);

% kontrolne točke Coonsove ploskve
[Bx,By,Bz] = coons(Bx,By,Bz)

% 2)
[u,v] = deal(linspace(0,1,50));
plotbezier2(Bx,By,Bz,u,v);
