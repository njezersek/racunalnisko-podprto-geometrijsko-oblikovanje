function d = deljenadif(r,b,i)
% Izračuna r-to deljeno diferenco b_i, kjer je b matrika ali stolpični vektor

if r == 0 % todo: dodaj d = 0, če je i izven meja
    d = b(i,:);
else
    d = deljenadif(r-1, b, i+1) - deljenadif(r-1, b, i);
end

end

