import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def pointbary(P, T):
# function u = pointbary(P,T)
# % Opis:
# % pointbary vrne baricentrične koordinate točke glede na dan trikotnik
# %
# % Definicija:
# % u = pointbary(P,T)
# %
# % Vhodna podatka:
# % P vrstica, ki predstavlja kartezične koordinate točke v ravnini,
# % T matrika s tremi vrsticami, v kateri vsaka vrstica predstavlja
# % kartezične koordinate oglišč trikotnika
# %
# % Izhodni podatek:
# % u vrstica dolžine 3, ki predstavlja baricentrične koordinate
# % točke, podane s P, glede na trikotnik, podan s T
# end

    A = np.array([[1, 1, 1], 
                  [T[0,0], T[1,0], T[2,0]],
                  [T[0,1], T[1,1], T[2,1]]])

    b = np.array([[1], [P[0]], [P[1]]])

    x = np.linalg.solve(A,b)
    # x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)

    return (x.T)[0]

def vectorbary(P, T):
    # enako kot pointbary, le da je b[0] = 0 (baje se poracuna uizi, da prides do tega)    
    A = np.array([[1, 1, 1], 
                  [T[0,0], T[1,0], T[2,0]],
                  [T[0,1], T[1,1], T[2,1]]])

    b = np.array([[0], [P[0]], [P[1]]])

    x = np.linalg.solve(A,b)
    # x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)

    return (x.T)[0]

def blossom3(B,U):
# function b = blossom3(B,U)
# % Opis:
# % blossom3 izračuna razcvet polinoma dveh spremenljivk
# %
# % Definicija:
# % b = blossom3(B,U)
# %
# % Vhodna podatka:
# % B matrika velikosti n+1 x n+1, ki predstavlja koeficiente
# % polinoma dveh spremenljivk stopnje n v Bezierjevi obliki
# % (element matrike na mestu (i,j), j <= n+2-i, določa
# % koeficient polinoma z indeksom (n+2-i-j, j-1, i-1)),
# % U matrika velikosti n x 3, v kateri vrstice predstavljavo
# % baricentrične koordinate točk glede na domenski trikotnik,
# % za katere izvajamo razcvet polinoma
# %
# % Izhodni podatek:
# % b vrednost razcveta polinoma, določenega z matriko B, v
# % točkah, določenih z matriko U
# end

# DeCasteljeujev postopek je to isto, le da namesto U[i] na vsakem koraku uporabimo enako trojko (u, v, w)

    A = np.copy(B)
    n = (U.shape)[0]
    
    for r in range(n):
        u = U[r, 0]
        v = U[r, 1]
        w = U[r, 2]
        
        for i in range(n - r):
            for j in range(n-r-i):
                A[i,j] = u * A[i,j] + v * A[i, j+1] + w * A[i+1, j]

    return A[0,0]

def decasteljau3(B, U):
    # enako kot blossom3, le da je U 1x3, uporabljas vedno iste u,v,w
    n = (B.shape)[0] - 1
    return blossom3(B, np.tile(U, (n, 1))) # tile ti multiplicira np array n-krat

def trimeshgrid(d):
    """
    Generates a triangulation of a triangle with vertices (0,0), (1,0), and (0,1).
    
    Args:
        d (int): Number of intervals dividing each side of the triangle.
    
    Returns:
        tuple: (TRI, U) where TRI is an array of triangles and U is an array of barycentric coordinates.
    """
    TRI = []
    U = []
    k = 0
    for i in range(d + 1):
        for j in range(d + 1 - i):
            U.append([d - i - j, j, i])
            if i < d and j < d - i:
                shift = d + 1 - i
                TRI.append([k, k + 1, k + shift])
                if j < d - i - 1:
                    TRI.append([k + 1, k + shift, k + 1 + shift])
            k += 1
    U = np.array(U) / d  # Normalize barycentric coordinates
    TRI = np.array(TRI)
    return TRI, U

def bezier3(Bx, By, Bz, U):
# function b = bezier3(Bx,By,Bz,U)
# % Opis:
# % bezier3 izračuna točke na trikotni Bezierjevi ploskvi
# %
# % Definicija:
# % b = bezier3(Bx,By,Bz,U)
# %
# % Vhodni podatki:
# % Bx, By, Bz matrike velikosti n+1 x n+1, ki predstavljajo
# % koordinate kontrolnih točk Bezierjeve krpe (element
# % posamezne matrike na mestu (i,j), j <= n+2-i, določa
# % koordinato kontrolne točke z indeksom
# % (n+2-i-j, j-1, i-1)),
# % U matrika, v kateri vrstice predstavljajo baricentrične
# % koordinate točk glede na domenski trikotnik, za katere
# % računamo točke na Bezierjevi krpi
# %
# % Izhodni podatek:
# % b matrika, v kateri vsaka vrstica predstavlja točko na
# % Bezierjevi krpi pri istoležnih parametrih iz matrike U
# end

    k = U.shape[0]

    b = np.zeros((k, 3))

    for i, u in enumerate(U):
        b[i] = decasteljau3(np.array([Bx, By, Bz]), u)
        # b[i, 1] = decasteljau3(By, u)
        # b[i, 2] = decasteljau3(Bz, u)

    return b

def plotbezier3(Bx, By, Bz, d):

    TRI, U = trimeshgrid(d)
    
    b = bezier3(Bx, By, Bz, U)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    ax.plot_trisurf(b[:, 0], b[:, 1], b[:, 2], triangles=TRI, alpha=0.8, cmap='viridis', edgecolor='none')

    # Plot the control mesh
    n = Bx.shape[0]
    for i in range(n):
        for j in range(n - i):
            if i < n - 1:
                ax.plot([Bx[i, j], Bx[i + 1, j]], [By[i, j], By[i + 1, j]], [Bz[i, j], Bz[i + 1, j]], 'k-')
            if j < n - i - 1:
                ax.plot([Bx[i, j], Bx[i, j + 1]], [By[i, j], By[i, j + 1]], [Bz[i, j], Bz[i, j + 1]], 'k-')
    ax.scatter(Bx.flatten(), By.flatten(), Bz.flatten(), color='r')
    plt.show()

if __name__ == '__main__':
    T = np.array([[0, 0], [5, 1], [3, 3]])
    p1 = pointbary([0, 0], T)
    print(p1)
    p2 = pointbary([1, 1], T)
    print(p2)
    p3 = pointbary([4, 2], T)
    print(p3)
    
    x = vectorbary([1, 0], T)
    y = vectorbary([0, 1], T)
    print(x, y)
    
    B = np.array([[2, 1, -1, 0],
                  [5, 3, -4, np.nan],
                  [0, 0, np.nan, np.nan],
                  [1, np.nan, np.nan, np.nan]])
    
    print(np.array([decasteljau3(B, p1), decasteljau3(B, p2), decasteljau3(B, p3)]))
    
    print(3*np.array([blossom3(B,np.array([x,p1,p1])), blossom3(B,np.array([x,p2,p2])), blossom3(B,np.array([x,p3,p3]))]))
    
    print(3*np.array([blossom3(B,np.array([y,p1,p1])), blossom3(B,np.array([y,p2,p2])), blossom3(B,np.array([y,p3,p3]))]))
    
    print(6*np.array([blossom3(B,np.array([x,x,p1])), blossom3(B,np.array([x,x,p2])), blossom3(B,np.array([x,x,p3]))]))
    
    print(6*np.array([blossom3(B,np.array([x,y,p1])), blossom3(B,np.array([x,y,p2])), blossom3(B,np.array([x,y,p3]))]))
    
    print(6*np.array([blossom3(B,np.array([y,y,p1])), blossom3(B,np.array([y,y,p2])), blossom3(B,np.array([y,y,p3]))]))
    
    # Bx = np.array([[4, 2, 5, 8], [0, 1, 7, 0], [0, 2, 6, 0], [0, 0, 0, 0]])
    # By = np.array([[5, 4, 3, -1], [0, 2, 3, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    # Bz = np.array([[3, 0, 5, 0], [-2, 4, 3, 0], [0, -2, 0, 0], [0, 0, 0, 0]])

    Bx = np.array([[0, 2, 6, 8], [1, 3, 7, 0], [2, 5, 0, 0], [4, 0, 0, 0]])
    By = np.array([[0, 1, 0, -1], [2, 2, 3, 0], [4, 3, 0, 0], [5, 0, 0, 0]])
    Bz = np.array([[-2, 1, -2, 0], [4, -2, 3, 0], [0, 5, 0, 0], [3, 0, 0, 0]])
    
    _, U = trimeshgrid(3)
    print(bezier3(Bx,By,Bz,U))

    # bezier3 neki zajebava, ampak je cudno da

    # Example usage
    # plotbezier3(Bx, By, Bz, d=50)