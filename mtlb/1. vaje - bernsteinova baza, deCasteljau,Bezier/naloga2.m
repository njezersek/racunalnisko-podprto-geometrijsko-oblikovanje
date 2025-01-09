% 2. NALOGA:

% 2.1)

% kontrolne točke Bezierjeve krivulje
B = [0 0; 1 2; 3 3; 4 -1; 0 2];

decasteljau(B(:,1),0)
decasteljau(B(:,2),0)
decasteljau(B(:,1),1)
decasteljau(B(:,2),1)
decasteljau(B(:,1),0.5)
decasteljau(B(:,2),0.5)

% 2.2)

B = [0 0; 1 2; 3 3; 4 -1; 0 2];
t = linspace(0,1,10);
bezier(B,t)

% 2.3)

B = [0 0; 1 2; 3 3; 4 -1; 0 2];
t = linspace(0,1);
figure
plotbezier(B,t);

% zrcaljenje krivulje čez y os
Bz = B*[-1 0; 0 1];
figure
plotbezier(Bz,t);

% rotacija krivulje za 60 stopinj
fi = pi/3;
Br = B*[cos(fi) -sin(fi); sin(fi) cos(fi)]';
figure
plotbezier(Br,t);
