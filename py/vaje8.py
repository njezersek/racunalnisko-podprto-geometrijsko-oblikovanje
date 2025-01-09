import numpy as np
from matplotlib import cm

from vaje2 import bezier
from vaje7 import bezier2

def lsqbezier2(m,n,P):
    # function B = lsqbezier2(m,n,P)
    # % Opis:
    # % lsqbezier2 vrne koeficiente tenzorskega polinoma, ki po metodi
    # % najmanjših kvadratov najbolje aproksimira dane podatke
    # %
    # % Definicija:
    # % B = lsqbezier2(m,n,P)
    # %
    # % Vhodni podatki:
    # % m,n parametra, ki določata stopnjo polinoma,
    # % P matrika podatkov, ki v vsaki vrstici vsebuje parametra z
    # % intervala [0,1] ter njima pripadajočo vrednost P[0,:] = (u0, v0, f0)
    # %
    # % Izhodni podatek:
    # % B matrika velikosti n+1 x m+1 s koeficienti polinoma, ki po
    # % metodi najmanjših kvadratov najbolje aproksimira podatke
    # end
    
    K = P.shape[0]
    
    A = np.zeros((K, (n+1)*(m+1)))
    
    for k in range(K):
        for i in range(m+1):
            for j in range(n+1):
                Bm = np.zeros((m+1,1))
                Bn = np.zeros((n+1,1))
                Bm[i] = 1
                Bn[j] = 1
                A[k, i*(n+1) + j] = bezier(Bm, [P[k,0]])[0] * bezier(Bn, [P[k,1]])[0]
    
    f = np.array(P[:,-1]).reshape(K, 1)
    # b = np.linalg.solve(A,f)
    x, residuals, rank, s = np.linalg.lstsq(A, f, rcond=None)
        
    return (x.reshape(m+1, n+1)).T

def lsqbezier2spline(M,N,m,n,P):
    # % Opis:
    # % lsqbezier2spline vrne kontrolne točke tenzorskih polinomov, ki
    # % določajo zvezno odvedljiv zlepek, ki v smislu metode najmanjših
    # % kvadratov najbolje aproksimira dane podatke
    # %
    # % Definicija:
    # % S = lsqbezier2spline(M,N,m,n,P)
    # %
    # % Vhodni podatki:
    # % M,N parametra, ki določata delitev domene parametrizacije,
    # % m,n parametra, ki določata stopnjo polinoma nad pravokotnikom
    # % delitve,
    # % P matrika podatkov, ki v vsaki vrstici vsebuje parametra z
    # % intervala [0,1] ter njima pripadajočo vrednost
    # %
    # % Izhodni podatek:
    # % S celica velikosti N x M, v kateri vsak element vsebuje
    # % MATRIKO (!) s koeficienti polinoma, ki določa zlepek nad
    # % pravokotnikom delitve domene
    # end
    
    print(P)
    
    S = np.zeros((N, M, (n+1), (m+1)))
    
    for I in range(M):
        for J in range(N): 
            trenutne_u = set(list(np.where(P[:,0] >= (I/M))[0])).intersection(set(list(np.where(P[:,0] <= ((I+1)/M))[0])))
            trenutne_v = set(list(np.where(P[:,1] >= (J/N))[0])).intersection(set(list(np.where(P[:,1] <= ((J+1)/N))[0])))
            presek = list(trenutne_u.intersection(trenutne_v))
            presek.sort()
            
            print(f"I: {I}, J: {J}")
            print(f"trenutne_u: {trenutne_u}")
            print(f"trenutne_v: {trenutne_v}")
            print(f"presek: {presek}")
                       
            if len(presek) > 0:
                P_sub = P[presek, :]
                print(f"P_sub: {P_sub}")
                S[J, I] = lsqbezier2(m, n, P_sub)
                print(f"S[{J}, {I}]: {S[J, I]}")
    
    return S


if __name__ == '__main__':
    
    # # V matlabu
    # # [X,Y,Z] = peaks(50);
    # # P = [(X(:)+3)/6 (Y(:)+3)/6 Z(:)];
    x = np.linspace(-3, 3, 50)
    y = np.linspace(-3, 3, 50)
    X, Y = np.meshgrid(x, y)
    Z = 3 * (1 - X)**2 * np.exp(-(X**2) - (Y+1)**2) - \
        10 * (X/5 - X**3 - Y**5) * np.exp(-X**2 - Y**2) - \
        1/3 * np.exp(-(X+1)**2 - Y**2)

    # Flatten and transform data
    P = np.column_stack(((X.ravel() + 3) / 6, (Y.ravel() + 3) / 6, Z.ravel()))
        
    # m = 6
    # n = 5
    # [Bx,By] = np.meshgrid(np.linspace(-3,3,m+1),np.linspace(-3,3,n+1))
    # Bz = lsqbezier2(m,n,P)
    
    # # print(Bz)
    
    # # # [u,v] = deal(linspace(0,1,7))
    # # # [~,~,bz] = bezier2(Bx,By,Bz,u,v)

    # u = np.linspace(0, 1, 7)
    # v = np.linspace(0, 1, 7)
    
    # _, _, bz = bezier2(Bx, By, Bz, u, v)
    # # print(bz)

    # u = np.linspace(0,1,50)
    # v = np.linspace(0,1,50)
    # _, _, bz = bezier2(Bx, By, Bz, u, v) 
    
    # print(np.linalg.norm(bz.reshape(-1) - P[:,-1],np.inf))
    
    # m = 9 
    # n = 10
    # [Bx,By] = np.meshgrid(np.linspace(-3,3,m+1),np.linspace(-3,3,n+1))
    # Bz = lsqbezier2(m,n,P)
    # [bx,by,bz] = bezier2(Bx,By,Bz,u,v)
    # print(np.linalg.norm(bz.reshape(-1) - P[:,-1],np.inf))

    # ZGORAJ VSE DELA! 
    # SPODAJ PA NE!

    # M = 5
    # N = 4
    # m = 3
    # n = 4
    # Sz = lsqbezier2spline(M, N, m, n, P)
    
    # print(Sz)