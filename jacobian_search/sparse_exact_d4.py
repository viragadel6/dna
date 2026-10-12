import sympy as sp, random, itertools

x,y=sp.symbols('x y')

monoms_d4 = [(2,0),(1,1),(0,2),(3,0),(2,1),(1,2),(0,3),(4,0),(3,1),(2,2),(1,3),(0,4)]

def monom_expr(i,j):
    return x**i * y**j

# Random sampling of supports size 4 for A and B
# For each support, build system and try to solve exactly via Groebner or solve

def try_support(SA, SB):
    # SA,SB list of (i,j)
    a_syms = sp.symbols(f'a0:{len(SA)}')
    b_syms = sp.symbols(f'b0:{len(SB)}')
    A = sum(a_syms[k]*monom_expr(*SA[k]) for k in range(len(SA)))
    B = sum(b_syms[k]*monom_expr(*SB[k]) for k in range(len(SB)))
    # Collision: A(1,0)= sum a where j=0
    A10 = sum(a_syms[k] for k in range(len(SA)) if SA[k][1]==0)
    B10 = sum(b_syms[k] for k in range(len(SB)) if SB[k][1]==0)
    eqs=[A10+1, B10]
    Ax=sp.diff(A,x)
    Ay=sp.diff(A,y)
    Bx=sp.diff(B,x)
    By=sp.diff(B,y)
    J=sp.expand(Ax+By+Ax*By-Ay*Bx)
    poly=sp.Poly(J, x,y)
    for coeff in poly.as_dict().values():
        eqs.append(coeff)
    # Try solve
    try:
        sol = sp.solve([e for e in eqs if e!=0], list(a_syms)+list(b_syms), dict=True)
        if sol:
            print(f"Found exact solution for SA={SA}, SB={SB}: {sol}")
            return True
    except Exception as e:
        pass
    # Try Groebner for inconsistency? Check if Groebner contains 1
    try:
        G = sp.groebner([e for e in eqs if e!=0], list(a_syms)+list(b_syms), order='lex')
        # If G contains 1, no solution
        # If G has no solution, skip
        # For small systems, Groebner may be heavy
        pass
    except:
        pass
    return False

# Random sampling
random.seed(0)
supports_d4 = list(itertools.combinations(monoms_d4, 4))
print(f"Total supports size4: {len(supports_d4)}")

for trial in range(100):
    SA = random.choice(supports_d4)
    SB = random.choice(supports_d4)
    # Ensure SA contains at least one pure x monomial for collision possibility
    if not any(j==0 for i,j in SA):
        continue
    if not any(j==0 for i,j in SB):
        continue
    print(f"Trying trial {trial}: SA={SA}, SB={SB}")
    found = try_support(SA, SB)
    if found:
        break

print("Done sparse exact d4")
