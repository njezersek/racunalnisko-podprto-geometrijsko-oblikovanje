addpath('..\1. vaje - bernsteinova baza, deCasteljau,Bezier', '..\2. vaje');

% 5. naloga:

% 1)
B = [0 0; 2 3; 4 2; 5 -1];
BS = beziersub(B,0.4);
[BS{1} BS{2}];

% 2)
BS1 = beziersub(B,0.5,1);
[BS1{1} BS1{2}];

BS2 = beziersub(B,0.5,2);
[BS2{1} BS2{2}; BS2{3} BS2{4}];

BS3 = beziersub(B,0.5,3);
[BS3{1} BS3{2}; BS3{3} BS3{4}; BS3{5} BS3{6}; BS3{7} BS3{8}];

% narišemo krivuljo pri različnih stopnjah subdivizije:
t = linspace(0,1);
subt = 0.5;
for k = 1:3
    hold on
    subplot(1,3,k);
    plotbezier(B,t);
    BS = beziersub(B,subt,k);
    for d = 1:length(BS)
        plotbezier(BS{d},t);
    end
end