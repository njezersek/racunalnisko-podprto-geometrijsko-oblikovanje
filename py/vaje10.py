import numpy as np
from scipy.spatial import Delaunay

from vaje9 import vectorbary

def f(V):
    x = V[0]
    y = V[1]
    return 3 * (1 - x)**2 * np.exp(-x**2 - (y + 1)**2) \
           - 10 * (x / 5 - x**3 - y**5) * np.exp(-x**2 - y**2) \
           - (1 / 3) * np.exp(-(x + 1)**2 - y**2)

def dxf(V):
    x = V[0]
    y = V[1]
    return (np.exp(-(x + 1)**2 - y**2) * (2 * x + 2)) / 3 \
           + 3 * np.exp(-(y + 1)**2 - x**2) * (2 * x - 2) \
           + np.exp(-x**2 - y**2) * (30 * x**2 - 2) \
           - 6 * x * np.exp(-(y + 1)**2 - x**2) * (x - 1)**2 \
           - 2 * x * np.exp(-x**2 - y**2) * (10 * x**3 - 2 * x + 10 * y**5)

def dyf(V):
    x = V[0]
    y = V[1]
    return (2 * y * np.exp(-(x + 1)**2 - y**2)) / 3 \
           + 50 * y**4 * np.exp(-x**2 - y**2) \
           - 3 * np.exp(-(y + 1)**2 - x**2) * (2 * y + 2) * (x - 1)**2 \
           - 2 * y * np.exp(-x**2 - y**2) * (10 * x**3 - 2 * x + 10 * y**5)

def Df(V):
    return np.array([dxf(V), dyf(V)])

def dxxf(V):
    x = V[0]
    y = V[1]
    return (2 * np.exp(-(x + 1)**2 - y**2)) / 3 \
           + 6 * np.exp(-(y + 1)**2 - x**2) \
           + 60 * x * np.exp(-x**2 - y**2) \
           - 6 * np.exp(-(y + 1)**2 - x**2) * (x - 1)**2 \
           - 2 * np.exp(-x**2 - y**2) * (10 * x**3 - 2 * x + 10 * y**5) \
           - (np.exp(-(x + 1)**2 - y**2) * (2 * x + 2)**2) / 3 \
           + 12 * x**2 * np.exp(-(y + 1)**2 - x**2) * (x - 1)**2 \
           + 4 * x**2 * np.exp(-x**2 - y**2) * (10 * x**3 - 2 * x + 10 * y**5) \
           - 12 * x * np.exp(-(y + 1)**2 - x**2) * (2 * x - 2) \
           - 4 * x * np.exp(-x**2 - y**2) * (30 * x**2 - 2)

def dxyf(V):
    x = V[0]
    y = V[1]
    return 4 * x * y * np.exp(-x**2 - y**2) * (10 * x**3 - 2 * x + 10 * y**5) \
           - (2 * y * np.exp(-(x + 1)**2 - y**2) * (2 * x + 2)) / 3 \
           - 2 * y * np.exp(-x**2 - y**2) * (30 * x**2 - 2) \
           - 100 * x * y**4 * np.exp(-x**2 - y**2) \
           - 3 * np.exp(-(y + 1)**2 - x**2) * (2 * x - 2) * (2 * y + 2) \
           + 6 * x * np.exp(-(y + 1)**2 - x**2) * (2 * y + 2) * (x - 1)**2

def dyyf(V):
    x = V[0]
    y = V[1]
    return (2 * np.exp(-(x + 1)**2 - y**2)) / 3 \
           - 6 * np.exp(-(y + 1)**2 - x**2) * (x - 1)**2 \
           - 2 * np.exp(-x**2 - y**2) * (10 * x**3 - 2 * x + 10 * y**5) \
           + 200 * y**3 * np.exp(-x**2 - y**2) \
           - 200 * y**5 * np.exp(-x**2 - y**2) \
           - (4 * y**2 * np.exp(-(x + 1)**2 - y**2)) / 3 \
           + 4 * y**2 * np.exp(-x**2 - y**2) * (10 * x**3 - 2 * x + 10 * y**5) \
           + 3 * np.exp(-(y + 1)**2 - x**2) * (2 * y + 2)**2 * (x - 1)**2

def Hf(V):
    return np.array([[dxxf(V), dxyf(V)], [dxyf(V), dyyf(V)]])

def argyris(T, f, Df, Hf):
    # % Opis:
    # % argyris izračuna Bezierjeve ordinate polinoma dveh spremenljivk
    # % stopnje 5, ki interpolira vrednosti, prve in druge odvode podane
    # % funkcije f v ogliščih trikotnika ter odvode funkcije f v središčih
    # % stranic trikotnika v smereh, pravokotnih na posamezno stranico
    # %
    # % Definicija:
    # % B = argyris(T,f,Df,Hf)
    # %
    # % Vhodni podatki:
    # % T tabela velikosti 3 x 2, v kateri vsaka vrstica predstavlja
    # % oglišče trikonika, nad katerim je definiran polinom,
    # % f funkcija, ki jo interpoliramo,
    # % Df gradient funkcije, ki jo interpoliramo,
    # % Hf Hessejeva matrika funkcije, ki jo interpoliramo
    # %
    # % Izhodni podatek:
    # % B matrika velikosti 6 x 6, ki predstavlja koeficiente polinoma
    # % dveh spremenljivk stopnje 5 v Bezierjevi obliki (element
    # % matrike na mestu (i,j), j <= 7-i, določa koeficient polinoma
    # % z indeksom (7-i-j, j-1, i-1)), ki interpolira funkcijo f
    # end
    
    # Initialize the coefficient matrix B
    B = np.zeros((6, 6))
    
        # Vertices of the triangle
    V1, V2, V3 = T[0], T[1], T[2]
    
    # 500
    B[0, 0] = f(V1)
    # 410
    B[0, 1] = f(V1) + (1/5) * np.dot(Df(V1), V2 - V1)
    # 401
    B[1, 0] = f(V1) + (1/5) * np.dot(Df(V1), V3 - V1)
    # 320
    B[0, 2] = f(V1) + (2/5) * np.dot(Df(V1), V2 - V1) + (1/20) * np.dot(V2 - V1, np.dot(Hf(V1), V2 - V1))
    # 311
    B[1, 1] = f(V1) + (1/5) * np.dot(Df(V1), V2 - V1) + (1/5) * np.dot(Df(V1), V3 - V1) + (1/20) * np.dot(V3 - V1, np.dot(Hf(V1), V2 - V1))
    # 302
    B[2, 0] = f(V1) + (2/5) * np.dot(Df(V1), V3 - V1) + (1/20) * np.dot(V3 - V1, np.dot(Hf(V1), V3 - V1))
    
    # 050
    B[0, 5] = f(V2)
    # 140
    B[0, 4] = f(V2) + (1/5) * np.dot(Df(V2), V1 - V2)
    # 041
    B[1, 4] = f(V2) + (1/5) * np.dot(Df(V2), V3 - V2)
    # 230
    B[0, 3] = f(V2) + (2/5) * np.dot(Df(V2), V1 - V2) + (1/20) * np.dot(V1 - V2, np.dot(Hf(V2), V1 - V2))
    # 131
    B[1, 3] = f(V2) + (1/5) * np.dot(Df(V2), V1 - V2) + (1/5) * np.dot(Df(V2), V3 - V2) + (1/20) * np.dot(V1 - V2, np.dot(Hf(V2), V3 - V2))
    # 032
    B[2, 3] = f(V2) + (2/5) * np.dot(Df(V2), V3 - V2) + (1/20) * np.dot(V3 - V2, np.dot(Hf(V2), V3 - V2))
    
    # 005
    B[5, 0] = f(V3)
    # 104
    B[4, 0] = f(V3) + (1/5) * np.dot(Df(V3), V1 - V3)
    # 014
    B[4, 1] = f(V3) + (1/5) * np.dot(Df(V3), V2 - V3)
    # 203
    B[3, 0] = f(V3) + (2/5) * np.dot(Df(V3), V1 - V3) + (1/20) * np.dot(V1 - V3, np.dot(Hf(V3), V1 - V3))
    # 113
    B[3, 1] = f(V3) + (1/5) * np.dot(Df(V3), V1 - V3) + (1/5) * np.dot(Df(V3), V2 - V3) + (1/20) * np.dot(V1 - V3, np.dot(Hf(V3), V2 - V3))
    # 023
    B[3, 2] = f(V3) + (2/5) * np.dot(Df(V3), V2 - V3) + (1/20) * np.dot(V2 - V3, np.dot(Hf(V3), V2 - V3))
    
    n12 = np.array([V2[1] - V1[1], V1[0] - V2[0]])  
    n12 = n12 / np.linalg.norm(n12)  
    theta = vectorbary(n12, T)
    
    V12 = 0.5*(V1 + V2)
    
    # 221
    B[1, 2] = (8/15) *(1/theta[2]) * np.dot(Df(V12), n12) - (1/6) * (B[1,0] + 4*B[1,1] + 4*B[1, 3] + B[1, 4]) \
              - (1/6) * (theta[0] / theta[2]) * (B[0, 0] + 4*B[0, 1] + 6*B[0, 2] + 4*B[0, 3] + B[0, 4]) \
              - (1/6) * (theta[1] / theta[2]) * (B[0, 5] + 4*B[0, 4] + 6*B[0, 3] + 4*B[0, 2] + B[0, 1])
    
    n23 = np.array([V3[1] - V2[1], V2[0] - V3[0]])  
    n23 = n23 / np.linalg.norm(n23)  
    theta = vectorbary(n23, T)
    
    V23 = 0.5*(V2 + V3)
    
    # 122
    B[2, 2] = (8/15) *(1/theta[0]) * np.dot(Df(V23), n23) - (1/6) * (B[0,4] + 4*B[1,3] + 4*B[3, 1] + B[4, 0]) \
              - (1/6) * (theta[1] / theta[0]) * (B[0, 5] + 4*B[1, 4] + 6*B[2, 3] + 4*B[3, 2] + B[4, 1]) \
              - (1/6) * (theta[2] / theta[0]) * (B[5, 0] + 4*B[4, 1] + 6*B[3, 2] + 4*B[2, 3] + B[1, 4])
    
    n13 = np.array([V3[1] - V1[1], V1[0] - V3[0]])  
    n13 = n13 / np.linalg.norm(n13) 
    theta = vectorbary(n13, T)
    
    V13 = 0.5*(V1 + V3)
    
    # 212
    B[2, 1] = (8/15) *(1/theta[1]) * np.dot(Df(V13), n13) - (1/6) * (B[0,1] + 4*B[1,1] + 4*B[3, 1] + B[4, 1]) \
              - (1/6) * (theta[0] / theta[1]) * (B[0, 0] + 4*B[1, 0] + 6*B[2, 0] + 4*B[3, 0] + B[4, 0]) \
              - (1/6) * (theta[2] / theta[1]) * (B[5, 0] + 4*B[4, 0] + 6*B[3, 0] + 4*B[2, 0] + B[1, 0])
    
    return B
    
def argyrisspline(tri,f,Df,Hf):
    # % Opis:
    # % argyrisspline izračuna Bezierjeve predstavitve polinomov, ki
    # % določajo Argyrisov zlepek nad triangulacijo
    # %
    # % Definicija:
    # % S = argyrisspline(tri,f,Df,Hf)
    # %
    # % Vhodni podatki:
    # % tri objekt razreda triangulation, ki določa triangulacijo
    # % domene, nad katero aproksimiramo funkcijo f,
    # % f funkcija, ki jo interpoliramo,
    # % Df gradient funkcije, ki jo interpoliramo,
    # % Hf Hessejeva matrika funkcije, ki jo interpoliramo
    # %
    # % Izhodni podatek:
    # % S celica z dolžino, ki ustreza številu trikotnikov v
    # % triangulaciji, v kateri vsak element vsebuje matriko
    # % velikosti 6 x 6, ki predstavlja koeficiente Argyrisovega
    # % polinoma v Bezierjevi obliki
    # end
    
    num_triangles = tri.simplices.shape[0]
    
    S = []
    
    for t in tri.simplices:
        T = tri.points[t]
        
        B = argyris(T, f, Df, Hf)
        
        S.append(B)
    
    return S
    
if __name__ == '__main__':
    T = np.array([[-1.7, -1], [-0.5, -1], [-0.5, 0]])

    Bz = np.round(argyris(T, f, Df, Hf),4)
    # print(Bz)
    
    # ARGYRIS DELA
        
    X = np.array([-3.0, -3.0, -3.0, -3.0, -1.7, -1.6, -1.5, -1.1, -0.9, -0.7, -0.5, -0.5, -0.5, -0.4, 0.0, 0.5, 0.5, 0.6, 0.7, 1.0, 1.2, 1.2, 1.6, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0])
    Y = np.array([-3.0, -2.0, 1.0, 3.0, -1.0, 0.5, 1.7, -3.0, -2.0, 1.0, -1.0, 0.0, 3.0, 1.6, -2.0, -1.0, 0.5, 1.2, -0.4, 2.2, -3.0, -1.6, 3.0, -1.0, 0.6, -3.0, -1.8, -0.4, 1.7, 3.0])
    tri = Delaunay(np.column_stack((X, Y)))
    S = argyrisspline(tri, f, Df, Hf)
    for i, B in enumerate(S):
        print(f"Triangle {i+1} Bezier coefficients matrix:\n{np.round(B, 4)}\n")
        
        # NJEGOVI S{1}, S{2}, S{3}, S{12} in S{15} 
        # SO MOJI B[25], B[?], B[?], B[15] in B[16]