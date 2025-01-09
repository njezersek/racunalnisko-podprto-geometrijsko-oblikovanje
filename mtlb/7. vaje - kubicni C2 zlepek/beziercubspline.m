function B = beziercubspline(u,D)
% Opis:
%   beziercubspline izračuna sestavljeno Bezierjevo krivuljo stopnje 3,
%   ki je dvakrat zvezno odvedljiva v stikih
%
% Definicija:
%   B = beziercubspline(u,D)
%
% Vhodna podatka:
%   u   seznam parametrov delitve dolžine m+1,
%   D   matrika, v kateri vsaka izmed m+3 vrstic predstavlja eno
%       kontrolno točko sestavljene krivulje
%
% Izhodni podatek:
%   B   celični seznam dolžine m, v kateri je vsak element matrika
%       s štirimi vrsticami, ki določajo kontrolne točke kosa
%       sestavljene krivulje

m = length(u) - 1;
[~,d] = size(D);
du = u(2:end) - u(1:end-1);
B = cell(1,m);

B{1} = zeros(4,d);
B{1}(1,:) = D(1,:); B{1}(2,:) = D(2,:);
B{1}(3,:) = (du(2)/(du(1)+du(2)))*D(2,:) + (du(1)/(du(1)+du(2)))*D(3,:);


for l=1:m-2
    B{l+1} = zeros(4,d);
    B{l+1}(2,:) = (du(l+1)+du(l+2))/(du(l)+du(l+1)+du(l+2))*D(l+2,:)+(du(l))/(du(l)+du(l+1)+du(l+2))*D(l+3,:);
    B{l+1}(3,:) = (du(l+2))/(du(l)+du(l+1)+du(l+2))*D(l+2,:)+(du(l)+du(l+1))/(du(l)+du(l+1)+du(l+2))*D(l+3,:);
end

B{m} = zeros(4,d);
B{m}(2,:) = du(m)/(du(m-1)+du(m))*D(m+1,:)+du(m-1)/(du(m-1)+du(m))*D(m+2,:);
B{m}(3,:) = D(m+2,:); B{m}(4,:) = D(m+3,:);

for l=1:m-1
    B{l}(4,:) = du(l+1)/(du(l)+du(l+1))*B{l}(3,:) + du(l)/(du(l)+du(l+1))*B{l+1}(2,:);
    B{l+1}(1,:) =  du(l+1)/(du(l)+du(l+1))*B{l}(3,:) + du(l)/(du(l)+du(l+1))*B{l+1}(2,:);
end

end