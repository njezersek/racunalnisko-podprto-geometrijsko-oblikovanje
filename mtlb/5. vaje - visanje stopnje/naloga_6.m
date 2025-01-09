addpath('..\1. vaje - bernsteinova baza, deCasteljau,Bezier', '..\2. vaje', '..\4. vaje - subdivizija');

% 6. naloga
B = [0 0; 2 6; 3 0; 5 4; 7 1; 5 5; 10 6];

bezierelv(B,1);
bezierelv(B,2);
bezierelv(B,3);

t = linspace(0,1);
k = [3 10 30];
for j = 1:3
    Br = bezierelv(B,k(j));
    subplot(1,3,j);
    plotbezier(B,t);
    plotbezier(Br,t);
end
