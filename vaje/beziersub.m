function BS = beziersub(B,t,k)
% Opis:
% beziersub izvede subdivizijo Bezierjeve krivulje
%
% Definicija:
% BS = beziersub(B,t)
%
% Vhodni podatki:
% B matrika kontrolnih točk Bezierjeve krivulje, v kateri
% vsaka vrstica predstavlja eno kontrolno točko,
% t parameter subdivizije Bezierjeve krivulje
% k parameter ki pove kolikokrat razpolovimo krivuljo
%
% Izhodni podatek:
% BS celica, ki vsebuje kontrolne točke dveh krivulj, ki jih
% dobimo s subdivizijo prvotne Bezierjeve krivulje

BS = {B};

for j = 1:k
    new_BS = {};
    for l = 1:length(BS)
        [n, m]  = size(BS{l});
        BL = zeros(n, m);
        BR = zeros(n, m);
        for i = 1:m
            D = decasteljau(BS{l}(:,i),t);
            BL(:,i) = D(1,:)';
            BR(:,i) = diag(fliplr(D));
        end
        
        new_BS = [new_BS, {BL, BR}];
    end
    BS = new_BS;
end

end
