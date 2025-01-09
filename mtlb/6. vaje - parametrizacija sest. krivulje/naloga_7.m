addpath('..\1. vaje - bernsteinova baza, deCasteljau,Bezier', '..\2. vaje', '..\4. vaje - subdivizija');

% 7. naloga

P = [-5 0; -2 -1; 0 3; 3 0; 7 -1];

% enakomerna parametrizacija:
ue = alphaparam(P,0)
% centripetalna parametrizacija:
uc = alphaparam(P,0.5)
% tetivna parametrizacija:
ut = alphaparam(P,1)

z = ones(1,5);
hold on
plot(ue,0*z,'ro-','MarkerEdgeColor','k','MarkerFaceColor','k');
plot(uc,z,'go-','MarkerEdgeColor','k','MarkerFaceColor','k');
plot(ut,2*z,'bo-','MarkerEdgeColor','k','MarkerFaceColor','k');
