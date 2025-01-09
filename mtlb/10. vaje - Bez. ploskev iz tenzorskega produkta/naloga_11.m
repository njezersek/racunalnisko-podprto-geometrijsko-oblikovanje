addpath('..\1. vaje - bernsteinova baza, deCasteljau,Bezier', '..\2. vaje', '..\4. vaje - subdivizija', '..\5. vaje - visanje stopnje','..\6. vaje - parametrizacija sest. krivulje','..\8. vaje - racionalne Bezierjeve krivulje');

% 11. naloga
% 1)

% kontrolne točke Bezierjeve ploskve
Bx = [0 3 5 6; 1 4 5 7; 0 3 5 6];
By = [-1 -2 -1 -2; 1 2 1 2; 5 5 6 5];
Bz = [-3 -2 0 -5; 3 6 3 2; -1 4 -2 8];

% izračun točk na ploskvi nad mrežo parametrov velikosti 7 x 7
[u,v] = deal(linspace(0,1,7));
[bx,by,bz] = bezier2(Bx,By,Bz,u,v);

% 2)

[u,v] = deal(linspace(0,1,50));
plotbezier2(Bx,By,Bz,u,v);
