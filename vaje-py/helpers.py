import numpy as np

# Naloga 2
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

	t = np.array(t)
	ta = np.array([t]).flatten() # poskrbimo, da je t 1D numpy array

	D = np.zeros((len(ta), dim, N, N)) * np.nan

	D[:,:,:,0] = B.T

	for r in range(1, N):
		D[:,:,:N-r,r] = ((1-ta)*D[:,:,:N-r,r-1].T + ta*D[:,:,1:N-r+1,r-1].T).T

	
	if dim == 1:
		D = D[:,0,:,:]
	
	if len(t.shape) == 0:
		D = D[0]

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