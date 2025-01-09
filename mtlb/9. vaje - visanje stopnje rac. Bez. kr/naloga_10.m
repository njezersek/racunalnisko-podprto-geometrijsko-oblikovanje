addpath('..\1. vaje - bernsteinova baza, deCasteljau,Bezier', '..\2. vaje', '..\4. vaje - subdivizija', '..\5. vaje - visanje stopnje','..\6. vaje - parametrizacija sest. krivulje','..\8. vaje - racionalne Bezierjeve krivulje');

% 10. naloga
% 1)
B = [1 0; 1 4; -3 2; -3 -2; 1 -4; 1 0];
w = [1; 1/5; 1/5; 1/5; 1/5; 1];

[B1,w1] = rbezierelv(B,w);
[B2,w2] = rbezierelv(B1,w1);

% 2) Preizkusite preprost postopek za aproksimacijo racionalne 
% Bézierjeve krivulje s polinomskimi, ki so določene s kontrolnimi 
% točkami racionalnih krivulj, dobljenih z višanjem stopnje prvotne 
% racionalne krivulje. Uspešnost postopka preverite na 
% parametrizaciji krožnice z racionalno krivuljo stopnje 5.

t = linspace(0,1);
b = rbezier(B,w,t);
N = 20;

clf
hold on

plot(b(:,1),b(:,2),'k','LineWidth',2);
C = jet(N);
for i = 1:N
    [B,w] = rbezierelv(B,w);
    b = bezier(B,t);
    plot(b(:,1),b(:,2),'Color',C(i,:));
end

hold off
axis equal
