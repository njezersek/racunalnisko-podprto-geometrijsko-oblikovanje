function S = lsqbezier2spline(M,N,m,n,P)
% Opis:
%   lsqbezier2spline vrne kontrolne točke tenzorskih polinomov, ki
%   določajo zvezno odvedljiv zlepek, ki v smislu metode najmanjših
%   kvadratov najbolje aproksimira dane podatke
%
% Definicija:
%   S = lsqbezier2spline(M,N,m,n,P)
%
% Vhodni podatki:
%   M,N     parametra, ki določata delitev domene parametrizacije,
%   m,n     parametra, ki določata stopnjo polinoma nad pravokotnikom
%           delitve,
%   P       matrika podatkov, ki v vsaki vrstici vsebuje parametra z
%           intervala [0,1] ter njima pripadajočo vrednost
%
% Izhodni podatek:
%   S       celica velikosti N x M, v kateri vsak element vsebuje
%           matriko s koeficienti polinoma, ki določa zlepek nad
%           pravokotnikom delitve domene

B = cell(N,M);

for I=1:M % stolpci
    for J=1:N % vrstice
        mask = (J/N>=P(:,2)).*(P(:,2)>=(J-1)/N).*(I/M>=P(:,1)).*(P(:,1)>=(I-1)/M);
        P_ = mask .* P;
        P_ = P_(any(P_,2),:);
        P_(:,1) = (P_(:,1)-(I-1)/M)*M; P_(:,2) = (P_(:,2)-(J-1)/N)*N;
        B{J,I} = lsqbezier2(m,n,P_);
    end
end

for I=1:M
    for J=1:N
        for i=0:m
            if J>1
                [B{J,I}(0+1,i+1),B{J-1,I}(n+1,i+1)] = deal(1/2.*B{J-1,I}(n-1+1,i+1) + 1/2.*B{J,I}(1+1,i+1));
            end
        end
        for j=0:n
            if I>1
                [B{J,I}(j+1,0+1),B{J,I-1}(j+1,m+1)] = deal(1/2.*B{J,I-1}(j+1,m-1+1) + 1/2.*B{J,I}(j+1,1+1));
            end
        end
        if I>1 && J>1
            [B{J,I}(0+1,0+1),B{J,I-1}(0+1,m+1),B{J-1,I-1}(n+1,m+1),B{J-1,I}(n+1,0+1)] = deal(1/4*B{J,I}(1+1,1+1) + 1/4*B{J,I-1}(1+1,m-1+1) + 1/4*B{J-1,I-1}(n-1+1,m-1+1) + 1/4*B{J-1,I}(n-1+1,1+1));
        end
    end
end

S = B;        

end