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


def test_rational_form(n):
	M = matrix.random(GF(257), n, n)
	while  M.minimal_polynomial().degree() != n:
		M = matrix.random(GF(257), n, n)
	print("M généré")
	start1 = walltime()
	F, U = rational_form_test(M, True)
	end1 = walltime()
	F2 = M.rational_form()
	end2 = walltime()
	F3 = rational_form_test(M, False)
	end3 = walltime()
	print("test d'égalité des formes: ", F2 == F == F3)
	print("test transformation ", U.inverse()*M*U == F)
	print("temps test avec transform: ", end1-start1, " secondes")
	print("temps de base: ", end2-end1, " secondes")
	print("temps test sans transform: ", end3-end2, " secondes")

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
	test, T = A.is_similar(B, transformation=True)
	end1 = walltime()
	test2, T2 = temps_nouveau(A, B)
	end2 = walltime()
	FA1 = A.rational_form()
	FB1 = B.rational_form()
	FA2 = rational_form_test(A, False)
	FB2 = rational_form_test(B, False)
	print("test d'égalité des matrices de passage: ", T == T2)
	print("test forme de Frobenius: ", FA1 == FA2 == FB1 == FB2)
	print("similarité 1: ", test)
	print("similarité 2: ", test2)
	print("test transformation ", T2.inverse()*B*T2 == A)
	print("temps de jordan: ", end1-start1, " secondes")
	print("temps cyclique: ", end2-end1, " secondes")

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