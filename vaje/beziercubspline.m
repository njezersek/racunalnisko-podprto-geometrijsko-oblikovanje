function B = beziercubspline(u,D)
% Opis:
% beziercubspline izračuna sestavljeno Bezierjevo krivuljo stopnje 3,
% ki je dvakrat zvezno odvedljiva v stikih
%
% Definicija:
% B = beziercubspline(u,D)
%
% Vhodna podatka:
% u seznam parametrov delitve dolžine m+1,
% D matrika, v kateri vsaka izmed m+3 vrstic predstavlja eno
% kontrolno točko sestavljene krivulje
%
% Izhodni podatek:
% B celični seznam dolžine m, v kateri je vsak element matrika
% s štirimi vrsticami, ki določajo kontrolne točke kosa
% sestavljene krivulje
m = length(u)-1;
d = size(D,2);
B = cell(m, 1);
for l = 1:m
    B{l} = zeros(4, d);
end


B{1}(1,:) = D(2-1,:); 
B{1}(2,:) = D(2+0,:); 
B{1}(3,:) = dd(u, 1+1)/(dd(u, 1+0) + dd(u, 1+1)) * D(2+0,:) + dd(u,1+0)/(dd(u, 1+0) + dd(u, 1+1)) * D(2+1,:);

for l = 1:m-2
    a = dd(u, 1+l-1) + dd(u, 1+l) + dd(u, 1+l+1);
    B{l+1}(2, :) = (dd(u, 1+l) + dd(u, 1+l+1))/a * D(2+l,:) + dd(u, 1+l-1) / a * D(2+l+1, :);
    B{l+1}(3, :) = dd(u, 1+l+1) / a * D(2+l,:) + (dd(u, 1+l-1) + dd(u, 1+l)) / a * D(2+l+1, :);
end

B{m}(2,:) = dd(u, 1+m-1)/(dd(u, 1+m-2)+dd(u, 1+m-1)) * D(2+m-1,:) + dd(u, 1+m-2) / (dd(u, 1+m-2)+dd(u, 1+m-1)) * D(2+m,:);
B{m}(3,:) = D(2+m,:);
B{m}(4,:) = D(2+m+1,:);

for l = 1:m-1
    a = dd(u, 1+l-1) + dd(u, 1+l);
    B{l}(4, :) = dd(u, 1+l) / a * B{l}(3,:) + dd(u, 1+l-1) / a * B{l+1}(2,:);
    B{l+1}(1, :) = B{l}(4, :);
end

end
