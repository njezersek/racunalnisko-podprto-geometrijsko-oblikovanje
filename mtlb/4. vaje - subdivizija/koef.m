function sez = koef(k)
%KOEF vrne seznam vrednosti, kolikokrat se mora določen interval
%razpoloviti pri računanju subdivizije

sez = zeros(1,2^k);
for i=1:k
    sez(2^(i-1)+1:2^i) = [sez(1:2^(i-1)-1) i];
end

sez = flip(sez);

end

