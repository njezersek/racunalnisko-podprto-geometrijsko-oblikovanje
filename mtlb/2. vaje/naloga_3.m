addpath('..\1. vaje - bernsteinova baza, deCasteljau,Bezier');

% 3. naloga:

B = [-2/3 -4/5; 1/3 1/5; 0 0; -1/3 1/5; 2/3 -4/5];
t = linspace(0,1,9);

% 3.1:
d1b = bezierder(B,1,t);
d2b = bezierder(B,2,t);
d3b = bezierder(B,3,t);
d4b = bezierder(B,4,t);

% 3.2:
[d1b,d1B] = bezierder(B,1,t);
[d2b,d2B] = bezierder(B,2,t);
[d3b,d3B] = bezierder(B,3,t);
[d4b,d4B] = bezierder(B,4,t);

% 3.3:
t = (0.46:0.01:0.54)';
d1b = bezierder(B,1,t);
tab = table();
tab{:,'t'} = t;
tab{:,'x'} = d1b(:,1);
tab{:,'y'} = d1b(:,2);
tab

t = linspace(0,1);
[d1b,d1B] = bezierder(B,1,t);
[d2b,d2B] = bezierder(B,2,t);
figure
plotbezier(B,t);
figure
plotbezier(d1B,t);
figure
plotbezier(d2B,t);
