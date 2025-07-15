def rational_form_test(M, transform):
	R = M.base_ring()
	n = M.nrows()
	poly_minimal = M.minimal_polynomial()
	if poly_minimal.degree() == n:
		if transform:
			v = vector(R, [1] + [0]*(n-1))
			iterates, _, _, _ = M._cyclic_subspace(v)
			T = matrix(iterates).transpose()
			return (T.inverse()*M*T, T)
		else:
			return companion_matrix(poly_minimal)
		
def rational_form_test_keller_gehrig(M, transform):
	if transform:
		return keller_gehrig(M, transformation=True, frobenius=True, poly=False)
	else:
		return keller_gehrig(M, frobenius=True, poly=False)

def keller_gehrig(A, transformation=False, poly=True, frobenius=False):
    R = A.base_ring()
    n = A.nrows()
    k = ceil(log(n, 2))
    v = vector(R, [1] + [0] * (n - 1))  # Vecteur unitaire e1

    M = A
    P = [v]
    for _ in range(k):
        P.extend([M * w for w in P])
        M = M*M

    # Construire la matrice U
    U = matrix(P).transpose()

    U = U[:, :n]
    if transformation and not poly and not frobenius:
        return U

    # Transformer A en forme de Frobenius simple
    U_inv = U.inverse()
    F = U_inv * A * U

    if frobenius and transformation and not poly:
        return (F, U)

    if frobenius and not transformation and not poly:
        return F

    # Extraire les coefficients du polynôme caractéristique
    char_poly = F.charpoly()
    coefficients = char_poly.coefficients(sparse=False)
    if frobenius and transformation and poly:
        return (coefficients, F, U_inv)

def krylov_iterative(A, v, transformation=False, poly=True, frobenius=False):
    R = A.base_ring()
    n = A.nrows()

    P = [v]
    for _ in range(n):
        P.append(A*P[-1])

    # Construire la matrice U
    U = matrix(P).transpose()

    U = U[:, :n]
    if transformation and not poly and not frobenius:
        return U

    # Transformer A en forme de Frobenius simple
    U_inv = U.inverse()
    F = U_inv * A * U

    if frobenius and transformation and not poly:
        return (F, U)

    if frobenius and not transformation and not poly:
        return F

    # Extraire les coefficients du polynôme caractéristique
    char_poly = F.charpoly()
    coefficients = char_poly.coefficients(sparse=False)
    if frobenius and transformation and poly:
        return (coefficients, F, U_inv)

def temps_nouveau_krylov(A, B):
	n = A.nrows()
	R = A.base_ring()
	v = vector(R, [1] + [0] * (n - 1))
	try:
		rat_form_A, SA = krylov_iterative(A, v, transformation=True, poly=False, frobenius=True)
		rat_form_B, SB = krylov_iterative(B, v, transformation=True, poly=False, frobenius=True)
		if rat_form_A == rat_form_B:
			return (True, SB * SA.inverse())
		else:
			return (False, None)
	except (ValueError, RuntimeError):
		pass

	try:
		ring = A.base_ring()
		closure = ring.algebraic_closure()
		A = A.change_ring(closure)
		B = B.change_ring(closure)
		rat_form_A, SA = rational_form_test(A, True)
		rat_form_B, SB = rational_form_test(B, True)
		if rat_form_A == rat_form_B:
			return (True, SB * SA.inverse())
		else:
			return (False, None)
	except (ValueError, RuntimeError, NotImplementedError):
		raise RuntimeError('unable to compute transformation for similar matrices')
	
def krylov_maxiaml_basis(A):
    F = A.base_ring()
    N = A.ncols() # En supposant que A est une matrice carrée n x n
    E = identity_matrix(N) # Base canonique
    
    U = matrix(F, N, N)
    k = 0

    for i in range(N):
        X = krylov_iterative(A, E.column(i), transformation=True, poly=False)
        Xi = X.transpose().rref().transpose()
        r = Xi.rank()
        U[:, k:r] = Xi
        print(Xi, end="\n\n")
        k += r 
    return U
	
def KGB(A):
    i = 0
    n = A.nrows()  # Supposant que A est une matrice carrée n x n
    V = [identity_matrix(A.base_ring(), n)] + [None] * (n-1)
    B = A

    while any(V[k].ncols() >= 2**i for k in range(i+1)):
        W = []
        k = 0
        for j in range(n):
            if  V[i][:, k:].ncols() < 2**i:
                Wj = V[i][:, k:]
                print("là")
            else:
                Wj = V[i].augment(B * V[i][:, k:])
                print("ici")
            W.append(Wj)
            print(Wj)
            k += V[i][k:].ncols()

        print(W)
        W_matrix = matrix(W)
        W_matrix_rref = W_matrix.transpose().rref().transpose()
        V.append(W_matrix_rref)
        B = B * B
        i += 1

    return V[i]

def KGBbis(A):
	R = M.base_ring()
	n = M.nrows()
	poly_minimal = M.minimal_polynomial()
	if poly_minimal.degree() == n:
		if transform:
			v = vector(R, [1] + [0]*(n-1))
			iterates, _, _, _ = M._cyclic_subspace(v)
			T = matrix(iterates).transpose()
			return (T.inverse()*M*T, T)
		else:
			return companion_matrix(poly_minimal)
		
def rational_form_test_keller_gehrig(M, transform):
	if transform:
		return keller_gehrig(M, transformation=True, frobenius=True, poly=False)
	else:
		return keller_gehrig(M, frobenius=True, poly=False)


def KGBbis(A):
    # Initialisation des matrices et vecteurs
    F = A.base_ring()
    N = A.ncols()
    U = matrix(F, N, N)  # Pour stocker A^2^i
    B = matrix(F, N, N)  # Pour stocker A^2^i
    V = matrix(F, N, N)  # Pour stocker A^2^i.U
    X = matrix(F, 2*N, N)  # Pour calculer la factorisation LSP

    P = list(range(N))  # Permutation de colonnes pour LQUP
    Q = list(range(2*N))  # Permutation de lignes pour LQUP

    d = [1] * N  # Dimensions de Vect(ei, Aei...)
    dv = d.copy()
    dold = d.copy()

    m = [[F(0) for _ in range(N)] for _ in range(N)]  # Opposé des coefficients des minpolys calculés
    i = 0
    l = 1
    k = N
    KeepOn = True

    # Initialisation de la matrice X avec (e1; e1A^t; e2; e2A^t;...;en;enA^t)
    for i in range(N):
        U.set_column(i, vector(F, [F.zero()] * N))
        U[i, i] = F.one()
        X.set_column(i, U.column(i))
        for j in range(N):
            X[N+i, j] = A[i*lda + j]
            B[i, j] = A[i*lda + j]
            V[i, j] = A[i*lda + j]

    # Factorisation LQUP
    LU, pivots = X.LU(decomposition=True, rowwise=False)

    k = update_d(F, d, KeepOn, l, N, X, Q, m)

    while KeepOn:  # Boucle principale
        # Mise à jour de U
        Uk = U
        cpt = 0
        for i in range(N):
            if d[i] < dold[i]:
                Ukp1new = Uk[d[i]:, :]
                Ukp1 = Uk[dold[i]:, :]
                Ukp1new[:] = Ukp1[:]
                dold[i] = d[i]
            cpt += d[i]

        # Insertion des blocs dupliqués
        Vk = V
        for i in range(k):
            Ukp1 = Uk[dold[i]:, :]
            newRowNb = d[i] - dold[i]
            if newRowNb > 0:
                Uk = U[-newRowNb:, :]
                Vk[:newRowNb, :] = Ukp1[:newRowNb, :]
            Uk = Uk[d[i]:, :]
            Vk = Vk[dv[i]:, :]

        # Mise à jour des dimensions
        l *= 2
        B = B^2

        # Calcul de V = U.B^t
        V = U * B.transpose()

        # Reconstruction de la matrice X
        X = block_matrix([[U[:cpt, :]], [V[:cpt, :]]])

        k = update_d(F, d, k, m)

        dv = d.copy()
        dold = d.copy()

        # Nouvelle factorisation LQUP
        LU, pivots = X.LU(decomposition=True, rowwise=False)

        # Mise à jour des dimensions
        k = update_d(F, d, KeepOn, l, N, X, Q, m)

    return charp


def temps_nouveau(A, B):
	try:
		rat_form_A, SA = rational_form_test(A, True)
		rat_form_B, SB = rational_form_test(B, True)
		if rat_form_A == rat_form_B:
			return (True, SB * SA.inverse())
		else:
			return (False, None)
	except (ValueError, RuntimeError):
		pass

	try:
		ring = A.base_ring()
		closure = ring.algebraic_closure()
		A = A.change_ring(closure)
		B = B.change_ring(closure)
		rat_form_A, SA = rational_form_test(A, True)
		rat_form_B, SB = rational_form_test(B, True)
		if rat_form_A == rat_form_B:
			return (True, SB * SA.inverse())
		else:
			return (False, None)
	except (ValueError, RuntimeError, NotImplementedError):
		raise RuntimeError('unable to compute transformation for similar matrices')

def temps_nouveau_keller(A, B):
	try:
		rat_form_A, SA = rational_form_test_keller_gehrig(A, True)
		rat_form_B, SB = rational_form_test_keller_gehrig(B, True)
		if rat_form_A == rat_form_B:
			print("erreur ici ?")
			return (True, SB * SA.inverse())
		else:
			return (False, None)
	except (ValueError, RuntimeError):
		pass

	try:
		ring = A.base_ring()
		closure = ring.algebraic_closure()
		A = A.change_ring(closure)
		B = B.change_ring(closure)
		rat_form_A, SA = rational_form_test_keller_gehrig(A, True)
		rat_form_B, SB = rational_form_test_keller_gehrig(B, True)
		if rat_form_A == rat_form_B:
			return (True, SB * SA.inverse())
		else:
			return (False, None)
	except (ValueError, RuntimeError, NotImplementedError):
		raise RuntimeError('unable to compute transformation for similar matrices')

def test_rational_form(n):
	M = matrix.random(GF(257), n, n)
	while  M.minimal_polynomial().degree() != n:
		M = matrix.random(GF(257), n, n)
	print("M généré")
	start1 = walltime()
	#F, U = rational_form_test(M, True)
	end1 = walltime()
	F2 = M.rational_form()
	end2 = walltime()
	F3 = keller_gehrig(M, frobenius=True, poly=False)
	end3 = walltime()
	F4, U2 = keller_gehrig(M, transformation=True, frobenius=True, poly=False)
	end4 = walltime()
	F5, U3 = keller_gehrig_alt(M, transformation=True, frobenius=True, poly=False)
	end5 = walltime()
	#print("test d'égalité des formes: ", F2 == F == F3 == F4)
	print("test transformation ", U2.inverse()*M*U2 == F4)
	print("test transformation 2: ", U3.inverse()*M*U3 == F5)
	print("temps test avec transform: ", end1-start1, " secondes")
	print("temps de base: ", end2-end1, " secondes")
	print("temps test keller-gehrig sans transform: ", end3-end2, " secondes")
	print("temps test keller-gehrig avec transform: ", end4-end3, " secondes")
	print("temps test keller-gehrig alt avec transform: ", end5-end4, " secondes")

def test_similar(n):
	A = matrix.random(GF(257), n, n)
	P = random_matrix(GF(257), n, n)
	while  A.minimal_polynomial().degree() != n:
		A = matrix.random(GF(257), n, n)
		while not P.is_invertible():
			P = random_matrix(GF(257), n, n)
	B = P.inverse()*A*P
	print("A et B généré")
	start1 = walltime()
	#test, T = A.is_similar(B, transformation=True)
	end1 = walltime()
	test2, T2 = temps_nouveau(A, B)
	end2 = walltime()
	test3, T3 = temps_nouveau_keller(A, B)
	end3 = walltime()
	#FA1 = A.rational_form()
	#FB1 = B.rational_form()
	FA2 = rational_form_test(A, False)
	FB2 = rational_form_test(B, False)
	FA3 = rational_form_test_keller_gehrig(A, False)
	FB3 = rational_form_test_keller_gehrig(B, False)
	print("test d'égalité des matrices de passage: ", T2 == T3)
	#print("test forme de Frobenius: ", FA1 == FA2 == FB1 == FB2 == FA3 == FB3)
	#print("similarité 1: ", test)
	print("similarité 2: ", test2)
	print("similarité 3: ", test3)
	print("test transformation ", T2.inverse()*B*T2 == A)
	print("test transformation keller-gehrig ", T3.inverse()*B*T3 == A)
	print("temps de jordan: ", end1-start1, " secondes")
	print("temps cyclique: ", end2-end1, " secondes")
	print("temps cyclique keller-gehrig: ", end3-end2, " secondes")

def test_max_similar():
	n = 100
	temps = 0
	while temps <= 1:
		n += 5
		A = matrix.random(GF(257), n, n)
		P = random_matrix(GF(257), n, n)
		while  A.minimal_polynomial().degree() != n:
			A = matrix.random(GF(257), n, n)
			while not P.is_invertible():
				P = random_matrix(GF(257), n, n)
		B = P.inverse()*A*P
		print("A et B généré pour n = ", n)
		start = walltime()
		_, T = A.is_similar(B, transformation=True)
		end = walltime()
		print("test transformation ", T.inverse()*B*T == A)
		temps = end-start
		print("temps pour n = ", n, " ", temps, " secondes")
		
		
def test_keller_gehrig(n):
	A = matrix.random(GF(257), n, n)
	P = random_matrix(GF(257), n, n)
	while  A.minimal_polynomial().degree() != n:
		A = matrix.random(GF(257), n, n)
		while not P.is_invertible():
			P = random_matrix(GF(257), n, n)
	B = P.inverse()*A*P
	print("A et B généré pour n = ", n)
	start = walltime()
	_, T = temps_nouveau_keller(A, B)
	end = walltime()
	print("test transformation ", T.inverse()*B*T == A)
	temps = end-start
	print("temps pour n = ", n, " ", temps, " secondes")
	
def test_KGB(n):
	A = matrix.random(GF(257), n, n)
	while  A.minimal_polynomial().degree() != n:
		A = matrix.random(GF(257), n, n)
	print("matrice généré")
	F = KGB(A)
	T = keller_gehrig_alt(A, transformation=True, poly=False)
	print(F.inverse()*A*F)
	print(T.inverse()*A*T)
	
def test_krylov_maximal(n):
    A = matrix.random(GF(257), n, n)
    print("matrice généré")
    U = krylov_maxiaml_basis(A)
    print(U.inverse()*A*U)
    # Initialisation des matrices et vecteurs
    F = A.base_ring()
    N = A.ncols()
    U = matrix(F, N, N)  # Pour stocker A^2^i
    B = matrix(F, N, N)  # Pour stocker A^2^i
    V = matrix(F, N, N)  # Pour stocker A^2^i.U
    X = matrix(F, 2*N, N)  # Pour calculer la factorisation LSP

    P = list(range(N))  # Permutation de colonnes pour LQUP
    Q = list(range(2*N))  # Permutation de lignes pour LQUP

    d = [1] * N  # Dimensions de Vect(ei, Aei...)
    dv = d.copy()
    dold = d.copy()

    m = [[F(0) for _ in range(N)] for _ in range(N)]  # Opposé des coefficients des minpolys calculés
    i = 0
    l = 1
    k = N
    KeepOn = True

    # Initialisation de la matrice X avec (e1; e1A^t; e2; e2A^t;...;en;enA^t)
    for i in range(N):
        U.set_column(i, vector(F, [F.zero()] * N))
        U[i, i] = F.one()
        X.set_column(i, U.column(i))
        for j in range(N):
            X[N+i, j] = A[i*lda + j]
            B[i, j] = A[i*lda + j]
            V[i, j] = A[i*lda + j]

    # Factorisation LQUP
    LU, pivots = X.LU(decomposition=True, rowwise=False)

    k = update_d(F, d, KeepOn, l, N, X, Q, m)

    while KeepOn:  # Boucle principale
        # Mise à jour de U
        Uk = U
        cpt = 0
        for i in range(N):
            if d[i] < dold[i]:
                Ukp1new = Uk[d[i]:, :]
                Ukp1 = Uk[dold[i]:, :]
                Ukp1new[:] = Ukp1[:]
                dold[i] = d[i]
            cpt += d[i]

        # Insertion des blocs dupliqués
        Vk = V
        for i in range(k):
            Ukp1 = Uk[dold[i]:, :]
            newRowNb = d[i] - dold[i]
            if newRowNb > 0:
                Uk = U[-newRowNb:, :]
                Vk[:newRowNb, :] = Ukp1[:newRowNb, :]
            Uk = Uk[d[i]:, :]
            Vk = Vk[dv[i]:, :]

        # Mise à jour des dimensions
        l *= 2
        B = B^2

        # Calcul de V = U.B^t
        V = U * B.transpose()

        # Reconstruction de la matrice X
        X = block_matrix([[U[:cpt, :]], [V[:cpt, :]]])

        k = update_d(F, d, k, m)

        dv = d.copy()
        dold = d.copy()

        # Nouvelle factorisation LQUP
        LU, pivots = X.LU(decomposition=True, rowwise=False)

        # Mise à jour des dimensions
        k = update_d(F, d, KeepOn, l, N, X, Q, m)

    return charp


def temps_nouveau(A, B):
	try:
		rat_form_A, SA = rational_form_test(A, True)
		rat_form_B, SB = rational_form_test(B, True)
		if rat_form_A == rat_form_B:
			return (True, SB * SA.inverse())
		else:
			return (False, None)
	except (ValueError, RuntimeError):
		pass

	try:
		ring = A.base_ring()
		closure = ring.algebraic_closure()
		A = A.change_ring(closure)
		B = B.change_ring(closure)
		rat_form_A, SA = rational_form_test(A, True)
		rat_form_B, SB = rational_form_test(B, True)
		if rat_form_A == rat_form_B:
			return (True, SB * SA.inverse())
		else:
			return (False, None)
	except (ValueError, RuntimeError, NotImplementedError):
		raise RuntimeError('unable to compute transformation for similar matrices')

def temps_nouveau_keller(A, B):
	try:
		rat_form_A, SA = rational_form_test_keller_gehrig(A, True)
		rat_form_B, SB = rational_form_test_keller_gehrig(B, True)
		if rat_form_A == rat_form_B:
			print("erreur ici ?")
			return (True, SB * SA.inverse())
		else:
			return (False, None)
	except (ValueError, RuntimeError):
		pass

	try:
		ring = A.base_ring()
		closure = ring.algebraic_closure()
		A = A.change_ring(closure)
		B = B.change_ring(closure)
		rat_form_A, SA = rational_form_test_keller_gehrig(A, True)
		rat_form_B, SB = rational_form_test_keller_gehrig(B, True)
		if rat_form_A == rat_form_B:
			return (True, SB * SA.inverse())
		else:
			return (False, None)
	except (ValueError, RuntimeError, NotImplementedError):
		raise RuntimeError('unable to compute transformation for similar matrices')

def test_rational_form(n):
	M = matrix.random(GF(257), n, n)
	while  M.minimal_polynomial().degree() != n:
		M = matrix.random(GF(257), n, n)
	print("M généré")
	start1 = walltime()
	#F, U = rational_form_test(M, True)
	end1 = walltime()
	F2 = M.rational_form()
	end2 = walltime()
	F3 = keller_gehrig(M, frobenius=True, poly=False)
	end3 = walltime()
	F4, U2 = keller_gehrig(M, transformation=True, frobenius=True, poly=False)
	end4 = walltime()
	F5, U3 = keller_gehrig_alt(M, transformation=True, frobenius=True, poly=False)
	end5 = walltime()
	#print("test d'égalité des formes: ", F2 == F == F3 == F4)
	print("test transformation ", U2.inverse()*M*U2 == F4)
	print("test transformation 2: ", U3.inverse()*M*U3 == F5)
	print("temps test avec transform: ", end1-start1, " secondes")
	print("temps de base: ", end2-end1, " secondes")
	print("temps test keller-gehrig sans transform: ", end3-end2, " secondes")
	print("temps test keller-gehrig avec transform: ", end4-end3, " secondes")
	print("temps test keller-gehrig alt avec transform: ", end5-end4, " secondes")

def test_similar(n):
	A = matrix.random(GF(257), n, n)
	P = random_matrix(GF(257), n, n)
	while  A.minimal_polynomial().degree() != n:
		A = matrix.random(GF(257), n, n)
		while not P.is_invertible():
			P = random_matrix(GF(257), n, n)
	B = P.inverse()*A*P
	print("A et B généré")
	start1 = walltime()
	#test, T = A.is_similar(B, transformation=True)
	end1 = walltime()
	test2, T2 = temps_nouveau(A, B)
	end2 = walltime()
	test3, T3 = temps_nouveau_keller(A, B)
	end3 = walltime()
	#FA1 = A.rational_form()
	#FB1 = B.rational_form()
	FA2 = rational_form_test(A, False)
	FB2 = rational_form_test(B, False)
	FA3 = rational_form_test_keller_gehrig(A, False)
	FB3 = rational_form_test_keller_gehrig(B, False)
	print("test d'égalité des matrices de passage: ", T2 == T3)
	#print("test forme de Frobenius: ", FA1 == FA2 == FB1 == FB2 == FA3 == FB3)
	#print("similarité 1: ", test)
	print("similarité 2: ", test2)
	print("similarité 3: ", test3)
	print("test transformation ", T2.inverse()*B*T2 == A)
	print("test transformation keller-gehrig ", T3.inverse()*B*T3 == A)
	print("temps de jordan: ", end1-start1, " secondes")
	print("temps cyclique: ", end2-end1, " secondes")
	print("temps cyclique keller-gehrig: ", end3-end2, " secondes")

def test_max_similar():
	n = 100
	temps = 0
	while temps <= 1:
		n += 5
		A = matrix.random(GF(257), n, n)
		P = random_matrix(GF(257), n, n)
		while  A.minimal_polynomial().degree() != n:
			A = matrix.random(GF(257), n, n)
			while not P.is_invertible():
				P = random_matrix(GF(257), n, n)
		B = P.inverse()*A*P
		print("A et B généré pour n = ", n)
		start = walltime()
		_, T = A.is_similar(B, transformation=True)
		end = walltime()
		print("test transformation ", T.inverse()*B*T == A)
		temps = end-start
		print("temps pour n = ", n, " ", temps, " secondes")
		
		
def test_keller_gehrig(n):
	A = matrix.random(GF(257), n, n)
	P = random_matrix(GF(257), n, n)
	while  A.minimal_polynomial().degree() != n:
		A = matrix.random(GF(257), n, n)
		while not P.is_invertible():
			P = random_matrix(GF(257), n, n)
	B = P.inverse()*A*P
	print("A et B généré pour n = ", n)
	start = walltime()
	_, T = temps_nouveau_keller(A, B)
	end = walltime()
	print("test transformation ", T.inverse()*B*T == A)
	temps = end-start
	print("temps pour n = ", n, " ", temps, " secondes")
	
def test_KGB(n):
	A = matrix.random(GF(257), n, n)
	while  A.minimal_polynomial().degree() != n:
		A = matrix.random(GF(257), n, n)
	print("matrice généré")
	F = KGB(A)
	T = keller_gehrig_alt(A, transformation=True, poly=False)
	print(F.inverse()*A*F)
	print(T.inverse()*A*T)
	
def test_krylov_maximal(n):
	A = matrix.random(GF(257), n, n)
	print("matrice généré")
	U = krylov_maxiaml_basis(A)
	print(U.inverse()*A*U)

def performance_krylov_rational(max):
	with open("test.txt", 'w') as f:
		for n in range(10, max, 10):
			A = matrix.random(GF(257), n, n)
			P = random_matrix(GF(257), n, n)
			while  A.minimal_polynomial().degree() != n:
				A = matrix.random(GF(257), n, n)
				while not P.is_invertible():
					P = random_matrix(GF(257), n, n)
			B = P.inverse()*A*P
			print("A et B généré pour n = ", n)
			start1 = walltime()
			#test, T = A.is_similar(B, transformation=True)
			end1 = walltime()
			test2, T2 = temps_nouveau(A, B)
			end2 = walltime()
			test3, T3 = temps_nouveau_krylov(A, B)
			end3 = walltime()
			f.write(f"Matrix size: {n}\n")
			#f.write(f"is_similar: {end1-start1}\n")
			f.write(f"Krylov Iterative Result: {end3-end2}\n")
			f.write(f"Rational Form Test Result: {end2-end1}\n")
			f.write("-" * 40 + "\n")
			print(f"Matrix size: {n}\n")
			#print(f"is_similar: {end1-start1}\n")
			print(f"Krylov Iterative Result: {end3-end2}\n")
			print(f"Rational Form Test Result: {end2-end1}\n")
			print(T2.inverse()*B*T2 == A and T3.inverse()*B*T2==A)
			print("-" * 40 + "\n")
		f.flush()
		
test_similar(2^10)