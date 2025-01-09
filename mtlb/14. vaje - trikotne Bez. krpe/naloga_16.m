addpath('..\13. vaje - polinomi dveh spremenljivk\');
% 16. naloga

% 1)
Bx = [0 2 6 8; 1 3 7 NaN; 2 5 NaN NaN; 4 NaN NaN NaN];
By = [0 1 0 -1; 2 2 3 NaN; 4 3 NaN NaN; 5 NaN NaN NaN];
Bz = [-2 1 -2 0; 4 -2 3 NaN; 0 5 NaN NaN; 3 NaN NaN NaN];

[~,U] = trimeshgrid(3);
bezier3(Bx,By,Bz,U)

% četrti parameter (50) določa vhodni podatek za metodo trimeshgrid
plotbezier3(Bx,By,Bz,50);
