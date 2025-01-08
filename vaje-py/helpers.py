import numpy as np
import matplotlib.pyplot as plt
import math

def matlab_to_numpy(matrix_str: str):
    """
    Converts a MATLAB-style matrix string to a NumPy array.
    
    Args:
        matrix_str (str): Matrix in MATLAB syntax (e.g., "[1 2; 3 4]").
    
    Returns:
        numpy.ndarray: Corresponding NumPy array.
    """
    # Remove MATLAB-style brackets
    matrix_str = matrix_str.strip("[]")
    
    # Split rows and columns
    rows = matrix_str.split(';')
    parsed_rows = []
    
    for row in rows:
        # Convert row string to a list of floats/integers
        parsed_rows.append([float(num) if '.' in num else int(num) for num in row.split()])
    
    # Convert the list of lists to a NumPy array
    return np.array(parsed_rows)

# Naloga 1. Zveza med Bernsteinovo in potenčno bazo

def power2bernstein(p: np.ndarray):
	"""
	Opis:
		power2bernstein pretvori polinom, predstavljen s koeficienti v
		potenčni bazi, v polinom, predstavljen v Bernsteinovi bazi

	Definicija:
		b = power2bernstein(p)

	Vhodni podatek:
		p	seznam koeficientov dolžine n+1, ki po vrsti pripadajo razvoju
			polinoma stopnje n v potenčni bazi od x^n do 1

	Izhodni podatek:
		b 	seznam koeficientov dolžine n+1, ki po vrsti pripadajo razvoju
			polinoma stopnje n v Bernsteinovi bazi od 0-tega do n-tega
			Bernsteinovega baznega polinoma
	"""	

	N = len(p)
	n = N-1
	b = np.zeros(N)

	for i in range(N):
		for j in range(i, N):
			b[j] = b[j] + math.comb(j,i)/math.comb(n,i) * p[n-i]

	return b


def bernstein2power(b: np.ndarray):
	"""
	Opis:
		bernstein2power pretvori polinom, predstavljen s koeficienti v
		bernsteinovi bazi, v polinom, predstavljen v potenčni bazi

	Definicija:
		p = bernstein2power(b)

	Vhodni podatek:
		b 	seznam koeficientov dolžine n+1, ki po vrsti pripadajo razvoju
			polinoma stopnje n v Bernsteinovi bazi od 0-tega do n-tega
			Bernsteinovega baznega polinoma

	Izhodni podatek:
		p 	seznam koeficientov dolžine n+1, ki po vrsti pripadajo razvoju
			polinoma stopnje n v potenčni bazi od x^n do 1
	"""

	N = len(b)
	n = N-1
	p = np.zeros(N)

	for i in range(N):
		for j in range(i, N):
			p[n-j] = p[n-j] + (-1)**(i+j) * math.comb(n, j) * math.comb(j, i) * b[i]
	
	return p

# Naloga 2. De Casteljaujev postopek

def decasteljau_basic(b, t):
	N = len(b)

	D = np.zeros((N,N)) * np.nan
	D[:,0] = b
	for r in range(1, N):
		D[:N-r,r] = (1-t)*D[:N-r,r-1] + t*D[1:N-r+1,r-1]

	return D

def decasteljau(B: np.ndarray, t: float | np.ndarray):
	"""
	Opis:
		decasteljau vrne shemo de Casteljaujevega postopka za dan seznam
		koordinat b pri danem parametru t
		
	Definicija:
		D = decasteljau(b,t)
		
	Vhodna podatka:
		b seznam koordinat kontrolnih točk Bezierjeve krivulje
		stopnje n,
		t parameter, pri katerem računamo koordinato Bezierjeve
		krivulje

	Izhodni podatek:
		D tabela velikosti n+1 x n+1, ki predstavlja de Casteljaujevo
		shemo za koordinate b pri parametru t (element na mestu
		(1,n+1) je koordinata Bezierjeve krivulje pri parametru t,
		elementi na mestih (i,j) za i > n-j+2 so NaN)
	"""

	N = len(B)
	dim = 1
	if len(B.shape) == 2:
		N, dim = B.shape
	if len(B.shape) > 2:
		raise ValueError("b mora biti 1D ali 2D numpy array")

	t = np.array([t]).flatten() # poskrbimo, da je t 1D numpy array

	D = np.zeros((len(t), dim, N, N)) * np.nan

	D[:,:,:,0] = B.T

	for r in range(1, N):
		D[:,:,:N-r,r] = ((1-t)*D[:,:,:N-r,r-1].T + t*D[:,:,1:N-r+1,r-1].T).T

	return D

def bezier(B, t):
	"""
	Opis:
		bezier vrne točke na Bezierjevi krivulji pri danih parametrih
	
	Definicija:
		b = bezier(B,t)
		
	Vhodna podatka:
		B matrika velikosti n+1 x d, ki predstavlja kontrolne točke
		Bezierjeve krivulje stopnje n v d-dimenzionalnem prostoru,
		t seznam parametrov dolžine k, pri katerih računamo vrednost
		Bezierjeve krivulje
	
	Izhodni podatek:
		b matrika velikosti k x d, kjer i-ta vrstica predstavlja točko
		na Bezierjevi krivulji pri parametru iz t na i-tem mestu
	"""

	return decasteljau(B, t)[:,:,0,-1]

def plotbezier(B: np.ndarray, t: np.ndarray, bezier_color: None | str = None, control_color = 'gray', controls=True):
	P = bezier(B, t)

	if P.shape[1] == 1:
		if controls: plt.plot(np.linspace(0,1,len(B)), B, 'o-', color=control_color, linewidth=1, alpha=0.6)
		plt.plot(t, P, '-', color=bezier_color, linewidth=2)
	if P.shape[1] == 2:
		if controls: plt.plot(B[:,0], B[:,1], 'o-', color=control_color, linewidth=1, alpha=0.6)
		plt.plot(P[:,0], P[:,1], '-', color=bezier_color, linewidth=2)
	if P.shape[1] == 3:
		ax = plt.figure().add_subplot(projection="3d")
		if controls: ax.plot(B[:,0], B[:,1], B[:,2], 'o-', color=control_color, linewidth=1, alpha=0.6)
		ax.plot(P[:,0], P[:,1], P[:,2], '-', color=bezier_color, linewidth=2)


# Naloga 3. Odvodi Bezierjeve krivulje

def deljenadif(r, b, i):
	"""
	Izračuna r-to deljeno diferenco b_i, kjer je b seznam (vektorjev)
	"""
	if r == 0:
		return b[i]
	else:
		return deljenadif(r-1, b, i+1) - deljenadif(r-1, b, i)
	
def bezierder(B, r, t):
	"""
		Opis:
			bezierder vrne točke na krivulji, ki predstavlja odvod dane
			Bezierjeve krivulje

		Definicija:
			db = bezierder(B,r,t)

		Vhodni podatki:
			B 	matrika kontrolnih točk Bezierjeve krivulje, v kateri vsaka
				vrstica predstavlja eno kontrolno točko,
			r 	stopnja odvoda, ki ga računamo,
			t 	seznam parameterov, pri katerih računamo odvod

		Izhodni podatek:
			db 	matrika, v kateri vsaka vrstica predstavlja točko r-tega
				odvoda pri istoležnem parametru iz seznama t

			dB 	matrika, v kateri vsaka vrstica predstavlja kontrolno točko
				r-tega odvoda
	"""
	B = np.array([[-2/3,-4/5],[1/3,1/5],[0,0],[-1/3,1/5],[2/3,-4/5]])
	N, d = B.shape
	n = N-1

	dB = np.zeros((N-r, d))
	for i in range(N-r):
		dB[i] = math.factorial(n)/math.factorial(n-r) * deljenadif(r, B, i)

	db = bezier(dB, t)

	return db, dB

# Naloga 4. Opis krožnice z Bezierjevimi krivuljami

def bezierarc(fi, m):
	"""
	Opis:
		bezierarc izračuna kontrole točke Bezierjeve krivulje, ki
		predstavlja interpolacijo krožnega loka po izbrani metodi

	Definicija:
		B = bezierarc(fi,m)

	Vhodna podatka:
		fi kot, ki določa krozni lok v območju (-fi,fi),
		m metoda interpolacije:
		1 = kvadratični G1 interpolant,
		2 = kubični C1 interpolant,
		3 = kubični G1 interpolant s C0 interpolacijo v srednji točki

	Izhodni podatek:
		B tabela velikosti 3 x 2 ali 4 x 2, v kateri vsaka vrstica
		predstavlja kontrolno točko Bezierjeve krivulje
	"""

	if m == 1: # kvadratični G1 interpolant
		B = np.zeros((3,2))
		B[0] = [np.cos(fi), -np.sin(fi)]
		B[1] = [1/np.cos(fi), 0]
		B[2] = [np.cos(fi), np.sin(fi)]
		return B
	
	if m == 2: # kubični C1 interpolant
		B = np.zeros((4,2))
		B[0] = [np.cos(fi), -np.sin(fi)]
		B[1] = [np.cos(fi)+1/3*np.sin(fi),-np.sin(fi)+1/3*np.cos(fi)]
		B[2] = [np.cos(fi)+1/3*np.sin(fi),np.sin(fi)-1/3*np.cos(fi)]
		B[3] = [np.cos(fi), np.sin(fi)]
		return B
	
	if m == 3: # kubični G1 interpolant s C0 interpolacijo v srednji točki
		B = np.zeros((4,2))
		B[0] = [np.cos(fi), -np.sin(fi)]
		B[1] = 1/3*np.array([4-np.cos(fi),4/np.tan(fi)-4/np.sin(fi)+np.sin(fi)])
		B[2] = 1/3*np.array([4-np.cos(fi),-4/np.tan(fi)+4/np.sin(fi)-np.sin(fi)])
		B[3] = [np.cos(fi), np.sin(fi)]
		return B

	raise ValueError("m mora biti 1, 2 ali 3")

# Naloga 5. Subdivizija Bezierjeve krivulje

def beziersub(B, t, k):
	"""
	Opis:
		beziersub izvede subdivizijo Bezierjeve krivulje
		
	Definicija:
		BS = beziersub(B,t)
		
	Vhodni podatki:
		B 	matrika kontrolnih točk Bezierjeve krivulje, v kateri
			vsaka vrstica predstavlja eno kontrolno točko,
		t 	parameter subdivizije Bezierjeve krivulje
		k 	parameter ki pove kolikokrat razpolovimo krivuljo
		
	Izhodni podatek:
		BS 	celica, ki vsebuje kontrolne točke dveh krivulj, ki jih
			dobimo s subdivizijo prvotne Bezierjeve krivulje
	"""

	BS = [B]

	for _ in range(k):
		new_BS = []
		for l in range(len(BS)):
			D = decasteljau(BS[l], t)
			BL = D[0,:,0,:].T
			BR = np.diagonal(np.fliplr(D[0, :]), axis1=1, axis2=2).T[::-1]
			new_BS += [BL, BR]
		BS = new_BS
	return BS


# Naloga 6. Višanje stopnje Bezierjeve krivulje

def bezierelf(B, k):
	pass

# Naloga 7
def alphaparam(P: np.ndarray, a: float):
	"""
	Opis:
		alphaparam sestavi alfa parametrizacijo oziroma delitev domene na
		podlagi točk

	Definicija:
		u = alphaparam(P,a)

	Vhodna podatka:
		P matrika z m+1 vrsticami, v kateri vsaka vrstica predstavlja
		eno točko,
		a parameter, ki določa alfa parametrizacijo

	Izhodni podatek:
		u seznam parametrov delitve, ki so določeni rekurzivno tako, da
		se trenutnemu parametru iz seznama u prišteje z a potencirana
		norma razlike zaporednih točk iz seznama P
	"""
	n, _ = P.shape
	u = np.zeros(n)
	u[0] = 0
	for i in range(1,n):
		u[i] = u[i-1] + np.linalg.norm(P[i] - P[i-1])**a
	return u


# Naloga 8
def beziercubspline(u, D):
	"""
	Opis:
		beziercubspline izračuna sestavljeno Bezierjevo krivuljo stopnje 3,
		ki je dvakrat zvezno odvedljiva v stikih

	Definicija:
		B = beziercubspline(u,D)

	Vhodna podatka:
		u seznam parametrov delitve dolžine m+1,
		D matrika, v kateri vsaka izmed m+3 vrstic predstavlja eno
		kontrolno točko sestavljene krivulje

	Izhodni podatek:
		B celični seznam dolžine m, v kateri je vsak element matrika
		s štirimi vrsticami, ki določajo kontrolne točke kosa
		sestavljene krivulje
	"""
	def d(u, i):
		if i < 0 or i+1 >= len(u):
			return 0
		return u[i+1] - u[i]

	D = np.roll(D, -1, axis=0) # popravimo indeksiranje v D

	m = len(u)-1
	dim = D.shape[1]

	B = [np.zeros((4, dim)) for _ in range(m)]


	B[0][0] = D[-1]
	B[0][1] = D[0]
	B[0][2] = d(u, 1) / (d(u, 0) + d(u,1)) * D[0] + d(u, 0) / (d(u, 0) + d(u,1)) * D[1]

	for l in range(1,m-1):
		a = d(u, l-1) + d(u, l) + d(u, l+1)
		B[l][1] = (d(u, l) + d(u, l+1)) / a * D[l] + d(u, l-1) / a * D[l+1]
		B[l][2] = d(u, l+1) / a * D[l] + (d(u, l-1) + d(u, l)) / a * D[l+1]

	B[m-1][1] = d(u, m-1) / (d(u, m-2) + d(u, m-1)) * D[m-1] + d(u, m-2) / (d(u, m-2) + d(u, m-1)) * D[m]
	B[m-1][2] = D[m]
	B[m-1][3] = D[m+1]

	for l in range(1,m):
		a = d(u, l-1) + d(u, l)
		B[l-1][3] = B[l][0] = d(u, l) / a * B[l-1][2] + d(u, l-1) / a * B[l][1]

	return B

# Naloga 9. Racionalne Bezierjeve krivulje in Farinove točke

def rbezier(b: np.ndarray, w: np.ndarray, t: np.ndarray):
	"""
	Opis:
		rbezier vrne točke na racionalni Bezierovi krivulji, izračunane z
		de Casteljaujevim postopkom za racionalne krivulje
		
	Definicija:
		b = rbezier(b,w,t)

	Vhod:
		b 	matrika velikosti n+1 x d, v kateri vsaka vrstica predstavlja
			eno kontrolno točko racionalne Bezierjeve krivulje stopnje n
		v 	prostoru dimenzije d,
		w 	seznam uteži racionalne Bezierjeve krivulje,
		t 	seznam parametrov dolžine N, za katere se računajo točke na
			racionalni Bezierjevi krivulji
		
	Izhod:
		 	matrika velikosti N x d, v kateri i-ta vrstica predstavlja
			točko na racionalni Bezierjevi krivulji pri i-tem parametru
			iz seznama t
	"""

	N, dim = b.shape

	B = np.zeros((len(t), dim, N, N)) # T, D, N, R
	W = np.zeros((len(t), 1, N, N)) # T, 1, N, R

	B[:,:,:,0] = b.T
	W[:,:,:,0] = w.T

	for r in range(1, N):
		W[:,:,:N-r,r] = ((1-t)*W[:,:,:N-r,r-1].T + t*W[:,:,1:N-r+1,r-1].T).T
		B[:,:,:N-r,r] = (1/W[:,:,:N-r,r].T * ((1-t) * W[:,:,:N-r,r-1].T * B[:,:,:N-r,r-1].T + t * W[:,:,1:N-r+1,r-1].T*B[:,:,1:N-r+1,r-1].T)).T

	return B[:, :, 0, -1]


def plotrbezier(B: np.ndarray, w: np.ndarray, t: np.ndarray, bezier_color: None | str = None, control_color = 'gray', controls=True):
	P = rbezier(B, w, t)

	if P.shape[1] == 1:
		if controls: plt.plot(np.linspace(0,1,len(B)), B, 'o-', color=control_color, linewidth=1, alpha=0.6)
		plt.plot(t, P, '-', color=bezier_color, linewidth=2)
		plt.gca().set_aspect('equal')
	if P.shape[1] == 2:
		if controls: 
			plt.plot(B[:,0], B[:,1], 'o-', color=control_color, linewidth=1, alpha=0.6)
			q = w[:-1]/(w[:-1]+w[1:])*B[:-1] + w[1:]/(w[:-1]+w[1:])*B[1:] # Farinove točke
			plt.plot(q[:,0], q[:,1], 'o', color='tab:blue', linewidth=1, alpha=0.6)
		plt.plot(P[:,0], P[:,1], '-', color=bezier_color, linewidth=2)

		plt.gca().set_aspect('equal')
	if P.shape[1] == 3:
		ax = plt.figure().add_subplot(projection="3d")
		if controls: ax.plot(B[:,0], B[:,1], B[:,2], 'o-', color=control_color, linewidth=1, alpha=0.6)
		ax.plot(P[:,0], P[:,1], P[:,2], '-', color=bezier_color, linewidth=2)