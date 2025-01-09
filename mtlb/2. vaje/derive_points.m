function dB = derive_points(B,r)
% dB = derive_points(B,r)
% Vrne dB-je, kjer je bezier(dB)=bezier'(B), za ' r-ti odvod in bezier(A) bezierovo krivuljo skozi tocke A.
    
    [n,d] = size(B); n = n-1;
    dB = zeros(n-r+1,d);
    coef = factorial(n)./factorial(n-r);

    for i=1:n-r+1
        dif = zeros(1,d);
        for j = 0:r
            dif = dif + nchoosek(r,j)*(-1)^(r-j).*B(i+j,:);
        end
        dB(i,:) = coef.*dif;
    end
end