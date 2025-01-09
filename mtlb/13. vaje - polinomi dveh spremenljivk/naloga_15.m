
% 15. naloga:

% 1)
% trikotnik
T = [0 0; 5 1; 3 3];
p1 = pointbary([0 0],T);
p2 = pointbary([1 1],T);
p3 = pointbary([4 2],T);

x = vectorbary([1 0],T);
y = vectorbary([0 1],T);

% 2)

% koeficienti polinoma v Bezierjevi obliki
B = [2 1 -1 0; 5 3 -4 NaN; 0 0 NaN NaN; 1 NaN NaN NaN];

% vrednosti p v točkah (0,0), (1,1), (4,2)
[decasteljau3(B,p1) decasteljau3(B,p2) decasteljau3(B,p3)];

% vrednosti Dxp v točkah (0,0), (1,1), (4,2)
3*[blossom3(B,[x;p1;p1]) blossom3(B,[x;p2;p2]) blossom3(B,[x;p3;p3])];

% vrednosti Dyp v točkah (0,0), (1,1), (4,2)
3*[blossom3(B,[y;p1;p1]) blossom3(B,[y;p2;p2]) blossom3(B,[y;p3;p3])];

% vrednosti Dxxp v točkah (0,0), (1,1), (4,2)
6*[blossom3(B,[x;x;p1]) blossom3(B,[x;x;p2]) blossom3(B,[x;x;p3])];

% vrednosti Dxyp v točkah (0,0), (1,1), (4,2)
6*[blossom3(B,[x;y;p1]) blossom3(B,[x;y;p2]) blossom3(B,[x;y;p3])];

% vrednosti Dyyp v točkah (0,0), (1,1), (4,2)
6*[blossom3(B,[y;y;p1]) blossom3(B,[y;y;p2]) blossom3(B,[y;y;p3])];