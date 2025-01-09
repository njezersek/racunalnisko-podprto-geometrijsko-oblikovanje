import numpy as np
import matplotlib.pyplot as plt
import math

def matlab(s: str):
	"""
	Converts a MATLAB-style matrix string to a NumPy array.
	
	Args:
		matrix_str (str): Matrix in MATLAB syntax (e.g., "[1 2; 3 4]").
	
	Returns:
		numpy.ndarray: Corresponding NumPy array.
	"""
	rows = s.split(";")
	res = "np.array(["+",".join([f'[{",".join(row.split())}]' for row in rows])+"])"
	a = eval(res)
	return a

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

def bezierelv(B, k):
	Be = B
	for j in range(k):
		n, m = Be.shape
		Bn = np.zeros((n+1, m))

		Bn[0] = Be[0]

		for i in range(1,n):
			Bn[i] = (1-i/n) * Be[i] + i/n * Be[i-1]
		
		Bn[n] = Be[n-1]

		Be = Bn
	return Be
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


def farin(B, w):
	"""Izračunaj farinove točke iz kontrolnih točk `B` in uteži `w`."""
	return w[:-1]/(w[:-1]+w[1:])*B[:-1] + w[1:]/(w[:-1]+w[1:])*B[1:]

def plotrbezier(B: np.ndarray, w: np.ndarray, t: np.ndarray, bezier_color: None | str = None, control_color = 'gray', controls=True):
	P = rbezier(B, w, t)

	if P.shape[1] == 1:
		if controls: plt.plot(np.linspace(0,1,len(B)), B, 'o-', color=control_color, linewidth=1, alpha=0.6)
		plt.plot(t, P, '-', color=bezier_color, linewidth=2)
		plt.gca().set_aspect('equal')
	if P.shape[1] == 2:
		if controls: 
			plt.plot(B[:,0], B[:,1], 'o-', color=control_color, linewidth=1, alpha=0.6)
			q = farin(B, w) # Farinove točke
			plt.plot(q[:,0], q[:,1], 'o', color='tab:blue', linewidth=1, alpha=0.6)
		plt.plot(P[:,0], P[:,1], '-', color=bezier_color, linewidth=2)

		plt.gca().set_aspect('equal')
	if P.shape[1] == 3:
		ax = plt.figure().add_subplot(projection="3d")
		if controls: ax.plot(B[:,0], B[:,1], B[:,2], 'o-', color=control_color, linewidth=1, alpha=0.6)
		ax.plot(P[:,0], P[:,1], P[:,2], '-', color=bezier_color, linewidth=2)


# Naloga 10. Višanje stopnje Bezierjeve krivulje

def rbezierelv(B, w, k=1):
	wBw = np.hstack([w*B, w])
	wBw_elv = bezierelv(wBw, k)

	w_elv = wBw_elv[:,-1, None]
	B_elv = wBw_elv[:,:-1] / w_elv

	return B_elv, w_elv 


# Naloga 11. Bezierjeve ploskve iz tenzorskega produkta.

def bezier2scaler(B, u, v):
	M, N = B.shape
	for _ in range(min(M, N)-1):
		M, N = B.shape
		B_new = np.zeros((M-1, N-1))
		for i in range(M-1):
			for j in range(N-1):
				B_new[i,j] = np.array([1-u, u]) @ B[i:i+2, j:j+2] @ np.array([1-v, v])
		B = B_new

	M, N = B.shape
	if N > M:
		for _ in range(N-M):
			M, N = B.shape
			B_new = np.zeros((M, N-1))
			for j in range(N-1):
				B_new[0,j] = (1-v) * B[0,j] + v * B[0,j+1]
			B = B_new
	if M > N:
		for _ in range(M-N):
			M, N = B.shape
			B_new = np.zeros((M-1, N))
			for i in range(M-1):
				B_new[i,0] = (1-u) * B[i,0] + u * B[i+1,0]
			B = B_new
	
	return float(B)

def bezier2(Bx, By, Bz, U, V):
	"""
	function [bx,by,bz] = bezier2(Bx,By,Bz,u,v)
	Opis:
		bezier2 vrne točke na Bezierjevi ploskvi iz tenzorskega produkta
	
	Definicija:
		[bx,by,bz] = bezier2(Bx,By,Bz,u,v)
	
	Vhodni podatki:
		Bx, By, Bz 	matrike velikosti n+1 x m+1, ki predstavljajo
					koordinate kontrolnih točk,
		u, v 		vrstici dolžine M in N, ki predstavljata parametre
					v smereh u in v
	
	Izhodni podatki:
		bx, by, bz matrike velikosti N x M, ki predstavljajo točke na
		Bezierjevi ploskvi: [bx(J,I) by(J,I) bz(J,I)] je
		točka pri parametrih u(I) in v(J)
	"""

	bx = np.zeros((len(U), len(V)))
	by = np.zeros((len(U), len(V)))
	bz = np.zeros((len(U), len(V)))
	for i, u in enumerate(U):
		for j, v in enumerate(V):
			bx[i, j] = bezier2scaler(Bx, u, v)
			by[i, j] = bezier2scaler(By, u, v)
			bz[i, j] = bezier2scaler(Bz, u, v)

	return bx, by, bz


def plotbezier2(Bx, By, Bz, U, V):
	bx, by, bz = bezier2(Bx, By, Bz, U, V)
	ax = plt.figure().add_subplot(projection='3d')
	ax.plot_surface(bx, by, bz, alpha=1, cmap="coolwarm") # type: ignore
	ax.plot_wireframe(Bx, By, Bz, alpha=0.5) # type: ignore
	ax.set_aspect('equal')



# Naloga 12. Coonsove ploskve

def coons1d(B):
	M, N = B.shape
	m = M-1
	n = N-1

	B1 = np.zeros((M, N))
	B2 = np.zeros((M, N))
	B3 = np.zeros((M, N))
	for i in range(M):
		for j in range(N):
			B1[i,j] = (1-j/n) * B[i, 0] + (j/n) * B[i, -1]
			B2[i,j] = (1-i/m) * B[0, j] + (i/m) * B[-1, j]
			B3[i,j] = (1-i/m)*(1-j/n)*B[0,0] + (1-i/m)*(j/n)*B[0,-1] + (i/m)*(1-j/n)*B[-1,0] + (i/m)*(j/n)*B[-1,-1]
			
	return B1 + B2 - B3
	

def coons(Bx, By, Bz):
	return coons1d(Bx), coons1d(By), coons1d(Bz)


# Naloga 13. Aprokicmacija z Bezierjevimi krivuljami po metodi najmanjših kvadratov

def peaks(n=50):
	X = np.linspace(-3, 3, n)
	Y = np.linspace(-3, 3, n)
	x, y = np.meshgrid(X, Y)
	z =  3*(1-x)**2 * np.exp(-(x**2) - (y+1)**2) \
		- 10*(x/5 - x**3 - y**5)*np.exp(-x**2-y**2) \
		- 1/3*np.exp(-(x+1)**2 - y**2)

	return x, y, z 

def lsqbezier2(m: int, n: int, P: np.ndarray):
	"""
	Opis:
		lsqbezier2 vrne koeficiente tenzorskega polinoma, ki po metodi
		najmanjših kvadratov najbolje aproksimira dane podatke

	Vhodni podatki:
		m,n parametra, ki določata stopnjo polinoma,
		P 	matrika podatkov, ki v vsaki vrstici vsebuje parametra z
			intervala [0,1] ter njima pripadajočo vrednost

	Izhodni podatek:
		B matrika velikosti n+1 x m+1 s koeficienti polinoma, ki po
		metodi najmanjših kvadratov najbolje aproksimira podatke	
	"""
	def bernstein(t, k, d):
		return (math.comb(d, k) * (t**k) * ((1 - t)**(d - k)))
	
	# Ekstrahiranje parametrov in vrednosti
	v = P[:, 0]
	u = P[:, 1]
	z = P[:, 2]
	
	# Dimenzije matrike A
	num_data = P.shape[0]
	num_coeffs = (m + 1) * (n + 1)
	
	# Ustvarimo matriko A
	A = np.zeros((num_data, num_coeffs))
	for i in range(n + 1):
		for j in range(m + 1):
			A[:, i * (m + 1) + j] = bernstein(u, i, n) * bernstein(v, j, m)
	
	# Rešitev problema najmanjših kvadratov za koeficiente
	B_vec, _, _, _ = np.linalg.lstsq(A, z, rcond=None)
	
	# Pretvorimo rezultat v matriko (n+1) x (m+1)
	B = B_vec.reshape((n + 1, m + 1))
	
	return B


# Naloga 14. Aproksimacija s sestavljenimi Bezierjevimi krivuljami

def lsqbezier2spline(M, N, m, n, P):
	I = np.arange(M+1)
	J = np.arange(N+1)

	U = I / M
	V = J / N

	B = dict()

	for i in range(1, M+1):
		for j in range(1, N+1):
			u = P[:, 0, None]
			v = P[:, 1, None]
			f = P[:, 2, None]

			mask = (U[i-1] <= u) & (u <= U[i]) & (V[j-1] <= v) & (v <= V[j])
			P_ = np.hstack(((u-U[i-1])/(U[i] - U[i-1]), (v-V[j-1])/(V[j]-V[j-1]), f))

			Bz = lsqbezier2(m, n, P_[mask.flatten()])
			B[(i, j)] = Bz

	# fix edges
	for II in range(1, M+1):
		for JJ in range(1, N+1):
			for ii in range(m+1):
				if JJ > 1:
					val = 1/2 * B[(II,JJ-1)][n-1, ii] + 1/2 * B[(II,JJ)][1, ii]
					B[(II, JJ)][0,ii] = B[(II,JJ-1)][n,ii] = val
			for jj in range(n+1):
				if II > 1:
					val = 1/2 * B[(II-1,JJ)][jj,m-1] + 1/2 * B[(II,JJ)][jj,1]
					B[(II, JJ)][jj,0] = B[(II-1,JJ)][jj,m] = val

			if II > 1 and JJ > 1:
				val = 1/4 * B[(II, JJ)][1,1] + 1/4 * B[(II-1,JJ)][1,m-1]  + 1/4 * B[(II-1,JJ-1)][n-1, m-1] + 1/4 * B[(II,JJ-1)][n-1,1]
				B[(II, JJ)][0,0] = val
				B[(II-1, JJ)][0,m] = val
				B[(II-1, JJ-1)][n,m] = val
				B[(II, JJ-1)][n,0] = val

	B_ = {(j, i): v for (i, j), v in B.items()} # obrni indekse v slovarju, da je isto kot v skripti

	return B_


### Naloga 15. Polinomi dveh spremenljivk

def pointbary(P: np.ndarray, T: np.ndarray):
	"""
	Opis:
		pointbary vrne baricentrične koordinate točke glede na dan trikotnik
	
	Vhodna podatka:
		P 	vrstica, ki predstavlja kartezične koordinate točke v ravnini,
		T	matrika s tremi vrsticami, v kateri vsaka vrstica predstavlja
			kartezične koordinate oglišč trikotnika
	
	Izhodni podatek:
		u vrstica dolžine 3, ki predstavlja baricentrične koordinate
		točke, podane s P, glede na trikotnik, podan s T
	"""
	A = np.array([[1, 1, 1], 
				[T[0,0], T[1,0], T[2,0]],
				[T[0,1], T[1,1], T[2,1]]])

	b = np.array([[1], [P[0]], [P[1]]])

	x = np.linalg.solve(A,b)

	return x.flatten()


def vectorbary(P, T):
	A = np.array([[1, 1, 1], 
				  [T[0,0], T[1,0], T[2,0]],
				  [T[0,1], T[1,1], T[2,1]]])

	b = np.array([[0], [P[0]], [P[1]]]) # tukaj je 0 namesto 1

	x = np.linalg.solve(A,b)

	return (x.T)[0]


def blossom3(B, U):
	"""
	Opis:
		blossom3 izračuna razcvet polinoma dveh spremenljivk

	Vhodna podatka:
		B 	matrika velikosti n+1 x n+1, ki predstavlja koeficiente
			polinoma dveh spremenljivk stopnje n v Bezierjevi obliki
			(element matrike na mestu (i,j), j <= n+2-i, določa
			koeficient polinoma z indeksom (n+2-i-j, j-1, i-1)),
		U	matrika velikosti n x 3, v kateri vrstice predstavljavo
			baricentrične koordinate točk glede na domenski trikotnik,
			za katere izvajamo razcvet polinoma

	Izhodni podatek:
		b 	vrednost razcveta polinoma, določenega z matriko B, v
			točkah, določenih z matriko U
	"""

	A = B.copy()	
	n = (U.shape)[0]
	
	for r in range(n):
		u = U[r, 0]
		v = U[r, 1]
		w = U[r, 2]
		
		for i in range(n - r):
			for j in range(n-r-i):
				A[i,j] = u * A[i,j] + v * A[i, j+1] + w * A[i+1, j]

	return A[0,0]


def decasteljau3(B, p):
	"""
		Opis:
			decasteljau3 izračuna vrednost polinoma, podanega s
			koeficienti B v točki p

		Vhodna podatka:
			B   matrika velikosti n+1 x n+1, ki predstavlja koeficiente
				polinoma dveh spremenljivk stopnje n v Bezierjevi obliki
				(element matrike na mestu (i,j), j <= n+2-i, določa
				koeficient polinoma z indeksom (n+2-i-j, j-1, i-1)),
			p   baricentrične koordinate točke, v kateri računamo vrednost polinoma

		Izhodni podatek:
			d	vrednost polinoma, določenega z B in p
	"""

	n = len(B)-1
	P = np.tile(p, (n, 1))
	return blossom3(B, P)

# Naloga 16. Trikotne Bezierjeve ploskve

def trimeshgrid(d):
	"""
	Opis:
		trimeshhgrid sestavi triangulacijo trikotnika z oglišči (0,0),
		(1,0), (0,1) in izračuna točke triangulacije v obliki
		baricentričnih koordinat

	Definicija:
	  	[TRI,U] = trimeshgrid(d)

	Vhodni podatek:
		d   število intervalov, na katere je razdeljena stranica
			trikotnika

	Izhodna podatka:
	  	TRI	matrika velikosti d^2 x 3, v kateri vsaka vrstica predstavlja
			en trikotnik enakomerne triangulacije, ki deli osnovni
			trikotnik,
		U	matrika velikosti (d+2)*(d+1)/2 x 3, v kateri vsaka vrstica
			predstavlja baricentrične koordinate ene točke triangulacije
	"""

	TRI = np.zeros((d**2, 3))
	U = np.zeros(((d+1)*(d+2)//2, 3))

	k = 0
	l = 0

	for i in range(d+1):
		for j in range(d-i+1):
			U[k] = np.array([d-i-j, j, i]) / d
			if i < d and j < d-i:
				shift = d+1-i
				TRI[l,:] = np.array([k, k+1, k+shift])
				l += 1
				if j < d-i-1:
					TRI[l,:] = np.array([k+1, k+shift, k+1+shift])
					l += 1
			k += 1

	return TRI, U

def bezier3(Bx, By, Bz, U):
	"""
	Opis:
		bezier3 izračuna točke na trikotni Bezierjevi ploskvi
	
	Vhodni podatki:
		Bx, By, Bz 	matrike velikosti n+1 x n+1, ki predstavljajo
			koordinate kontrolnih točk Bezierjeve krpe (element
			posamezne matrike na mestu (i,j), j <= n+2-i, določa
			koordinato kontrolne točke z indeksom
			(n+2-i-j, j-1, i-1)),
		U 	matrika, v kateri vrstice predstavljajo baricentrične
			koordinate točk glede na domenski trikotnik, za katere
			računamo točke na Bezierjevi krpi
	
	Izhodni podatek:
		b 	matrika, v kateri vsaka vrstica predstavlja točko na
			Bezierjevi krpi pri istoležnih parametrih iz matrike U

	"""
	
	k = len(U)
	b = np.zeros((k, 3))

	for i, u in enumerate(U):
		b[i, 0] = decasteljau3(Bx, u)
		b[i, 1] = decasteljau3(By, u)
		b[i, 2] = decasteljau3(Bz, u)
	
	return b


# Naloga 17. Argyrisova interpolacijska shema

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
    B = np.full((6, 6), np.nan)
    
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


# Naloga 18. Argyrisov zlepek

def argyrisspline(simplices: np.ndarray, points: np.ndarray, f, Df, Hf):
	"""
	Opis:
		argyrisspline izračuna Bezierjeve predstavitve polinomov, ki
		določajo Argyrisov zlepek nad triangulacijo

	Vhodni podatki:
		tri objekt razreda triangulation, ki določa triangulacijo
			domene, nad katero aproksimiramo funkcijo f,
		f 	funkcija, ki jo interpoliramo,
		Df	gradient funkcije, ki jo interpoliramo,
		Hf 	Hessejeva matrika funkcije, ki jo interpoliramo

	Izhodni podatek:
		S celica z dolžino, ki ustreza številu trikotnikov v
		triangulaciji, v kateri vsak element vsebuje matriko
		velikosti 6 x 6, ki predstavlja koeficiente Argyrisovega
		polinoma v Bezierjevi obliki
	"""

	S = []
	
	for t in simplices:
		T = points[t]
		
		B = argyris(T, f, Df, Hf)
		
		S.append(B)
	
	return S