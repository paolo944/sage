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

def keller_gehrig_alt(A, transformation=False, poly=True, frobenius=False):
    R = A.base_ring()
    n = A.nrows()
    v = vector(R, [1] + [0] * (n - 1))  # Vecteur unitaire e1

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
	
def KGB(A):
    i = 0
    n = A.nrows()  # Assuming A is a square matrix, n x n
    V = [identity_matrix(A.base_ring(), n)]  # V0 is the identity matrix In
    B = A

    while any(V[k].ncols() >= 2**i for k in range(n)):
        W = []
        for j in range(n):
            if V[i].column(j).ncols() < 2**i:
                Wj = V[i].column(j)
            else:
                Wj = V[i].augment(B * V[i].column(j))
            W.append(Wj)

        W_matrix = matrix(W)
        V.append(W_matrix.column_space().basis_matrix())
        B = B * B
        i += 1

    return V[i]

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
	
def test_keller_gehrig_non_cyclic(n):
	A = matrix.random(GF(2), n, n)
	while  A.minimal_polynomial().degree() == n:
		A = matrix.random(GF(2), n, n)
	F = KGB(A)
	print(F)