function [db, dB] = bezierder(B,r,t)
% Opis:
%   bezierder vrne točke na krivulji, ki predstavlja odvod dane
%   Bezierjeve krivulje
%
% Definicija:
%   [db, dB] = bezierder(B,r,t)
%
% Vhodni podatki:
%   B   matrika kontrolnih točk Bezierjeve krivulje, v kateri vsaka
%       vrstica predstavlja eno kontrolno točko,
%   r   stopnja odvoda, ki ga računamo,
%   t   seznam parameterov, pri katerih računamo odvod
%
% Izhodni podatek:
%   db  matrika, v kateri vsaka vrstica predstavlja točko r-tega
%       odvoda pri istoležnem parametru iz seznama t
%   dB  matrika kontrolnih točk Bézierjeve krivulje stopnje n − r, ki
%       predstavlja r-ti odvod podane Bézierjeve krivulje

[~,d] = size(B);

if nargout == 1
    
    n = length(B)-1;
    a = repmat([1 -1],1,r+1);
    a = a(1:r+1);
    koef = diag(fliplr(pascal(r+1)))'.*a;
    k = length(t);
    const = factorial(n)/factorial(n-r);
    db = zeros(k,d);
    
    if d == 2
        for i=1:k
            D1 = decasteljau(B(:,1),t(i));
            D2 = decasteljau(B(:,2),t(i));
            dB1 = const.*dot(flip(D1(1:r+1,n-r+1)),koef);
            dB2 = const.*dot(flip(D2(1:r+1,n-r+1)),koef);
            db(i,:) = [dB1 dB2];
        end
    elseif d == 3
        for i=1:k
            D1 = decasteljau(B(:,1),t(i));
            D2 = decasteljau(B(:,2),t(i));
            D3 = decasteljau(B(:,3),t(i));
            dB1 = const.*dot(flip(D1(1:r+1,n-r+1)),koef);
            dB2 = const.*dot(flip(D2(1:r+1,n-r+1)),koef);
            dB3 = const.*dot(flip(D3(1:r+1,n-r+1)),koef);
            db(i,:) = [dB1 dB2 dB3];
        end
    end

elseif nargout == 2

    dB = derive_points(B,r);
    k = length(t);
    db = zeros(k,d);
    if d == 2
        for i = 1:k
            D1 = decasteljau(dB(:,1),t(i));
            D2 = decasteljau(dB(:,2),t(i));
            db(i,:) = [D1(1,end),D2(1,end)];
        end
    elseif d == 3
        for i = 1:k
            D1 = decasteljau(dB(:,1),t(i));
            D2 = decasteljau(dB(:,2),t(i));
            D3 = decasteljau(dB(:,3),t(i));
            db(i,:) = [D1(1,end),D2(1,end),D3(1,end)];
        end
    end

end

end