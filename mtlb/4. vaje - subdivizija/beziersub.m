function BS = beziersub(B,t,k)
% Opis:
%   beziersub izvede subdivizijo Bezierjeve krivulje
%
% Definicija:
%   BS = beziersub(B,t)
%
% Vhodni podatki:
%   B   matrika kontrolnih točk Bezierjeve krivulje, v kateri
%       vsaka vrstica predstavlja eno kontrolno točko,
%   t   parameter subdivizije Bezierjeve krivulje
%   k   število zaporednih izvedenih subdivizij
%
% Izhodni podatek:
%   BS  celica, ki vsebuje kontrolne točke 2^k krivulj, ki jih
%       dobimo s k koraki subdivizije prvotne Bezierjeve krivulje

if nargin == 2
    k = 1;
end

BS = cell(1,2^k);
BS{1} = B;
sez = koef(k);

[~,d] = size(B);

while isempty(find(sez,1)) == 0
    i = find(sez,1);
    k = sez(i);
    B = BS{i};

    bx = B(:,1); by = B(:,2);
    vrednostix = decasteljau(bx,t);
    vrednostiy = decasteljau(by,t);
    kontrolne1x = vrednostix(1,:)';
    kontrolne2x = flip(diag(flip(vrednostix)));
    kontrolne1y = vrednostiy(1,:)';
    kontrolne2y = flip(diag(flip(vrednostiy)));
    kontrolne1 = [kontrolne1x kontrolne1y];
    kontrolne2 = [kontrolne2x kontrolne2y];

    BS{i} = kontrolne1; BS{i + 2^(k-1)} = kontrolne2;
    sez(i) = sez(i)-1;
end

if d == 3
    BS = cell(1,2^k);
    BS{1} = B;
    sez = koef(k);
    while isempty(find(sez,1)) == 0
        i = find(sez,1);
        k = sez(i);
        B = BS{i};
    
        bx = B(:,1); by = B(:,2); bz = B(:,3);
        vrednostix = decasteljau(bx,t);
        vrednostiy = decasteljau(by,t);
        vrednostiz = decasteljau(bz,t);
        kontrolne1x = vrednostix(1,:)';
        kontrolne2x = flip(diag(flip(vrednostix)));
        kontrolne1y = vrednostiy(1,:)';
        kontrolne2y = flip(diag(flip(vrednostiy)));
        kontrolne1z = vrednostiz(1,:)';
        kontrolne2z = flip(diag(flip(vrednostiz)));
        kontrolne1 = [kontrolne1x kontrolne1y kontrolne1z];
        kontrolne2 = [kontrolne2x kontrolne2y kontrolne2z];
    
        BS{i} = kontrolne1; BS{i + 2^(k-1)} = kontrolne2;
        sez(i) = sez(i)-1;
    end
end


end
