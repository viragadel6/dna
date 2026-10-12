import sympy as sp
import itertools

x,y=sp.symbols('x y')

# Monomials for d=3: degree 2..3
monoms = [(2,0),(1,1),(0,2),(3,0),(2,1),(1,2),(0,3)]
# Map to sympy monomial x^i y^j
def monom_to_expr(i,j):
    return x**i * y**j

# Enumerate supports of size 2 for A and B
# For each support, we need coefficients
# Let's try size 2 for each

supports = list(itertools.combinations(monoms, 2))
print(f"Number supports size2: {len(supports)}")

# For each pair of supports (SA,SB), build system and try solve
# We'll attempt for first 100 combos to save time

count=0
solutions_found=[]

for SA in supports:
    for SB in supports:
        count+=1
        if count>200:  # limit for now
            break
        # Build symbols
        # A = sum a_k * monom
        # B = sum b_k * monom
        a_syms = sp.symbols(f'a0:{len(SA)}')
        b_syms = sp.symbols(f'b0:{len(SB)}')
        A = sum(a_syms[i]*monom_to_expr(*SA[i]) for i in range(len(SA)))
        B = sum(b_syms[i]*monom_to_expr(*SB[i]) for i in range(len(SB)))
        # Check collision conditions: A(1,0) = sum a_i where j=0, i>=2
        A10 = sum(a_syms[i] for i in range(len(SA)) if SA[i][1]==0)
        B10 = sum(b_syms[i] for i in range(len(SB)) if SB[i][1]==0)
        # Equations: A10+1=0, B10=0
        eqs=[]
        eqs.append(A10+1)
        eqs.append(B10)
        # Jacobian condition
        Ax=sp.diff(A,x)
        Ay=sp.diff(A,y)
        Bx=sp.diff(B,x)
        By=sp.diff(B,y)
        J=sp.expand(Ax+By+Ax*By-Ay*Bx)
        poly=sp.Poly(J, x,y)
        for coeff in poly.as_dict().values():
            eqs.append(coeff)
        # Now we have system of eqs (len = 2 + up to 14) in 4 unknowns
        # Try to solve via Groebner or solve
        # Use sympy solve
        try:
            # Filter out trivial zero eqs
            eqs_nonzero = [e for e in eqs if e!=0]
            # If more eqs than vars, still try solve
            sol = sp.solve(eqs_nonzero, list(a_syms)+list(b_syms), dict=True)
            if sol:
                print(f"Found solution for SA={SA}, SB={SB}: {sol}")
                solutions_found.append((SA,SB,sol))
        except Exception as e:
            # print(f"Solve failed for {SA},{SB}: {e}")
            continue
    if count>200:
        break

print(f"Checked {count} combos, found {len(solutions_found)} solutions")

# Now try size 3 supports for d=3, random sampling

print("\n--- Random sampling size3 supports d=3 ---")
import random
random.seed(0)
monoms_d3 = monoms
supports3 = list(itertools.combinations(monoms_d3, 3))
print(f"Supports size3 count: {len(supports3)}")

for _ in range(100):
    SA = random.choice(supports3)
    SB = random.choice(supports3)
    a_syms = sp.symbols(f'a0:{len(SA)}')
    b_syms = sp.symbols(f'b0:{len(SB)}')
    A = sum(a_syms[i]*monom_to_expr(*SA[i]) for i in range(len(SA)))
    B = sum(b_syms[i]*monom_to_expr(*SB[i]) for i in range(len(SB)))
    A10 = sum(a_syms[i] for i in range(len(SA)) if SA[i][1]==0)
    B10 = sum(b_syms[i] for i in range(len(SB)) if SB[i][1]==0)
    eqs=[A10+1, B10]
    Ax=sp.diff(A,x)
    Ay=sp.diff(A,y)
    Bx=sp.diff(B,x)
    By=sp.diff(B,y)
    J=sp.expand(Ax+By+Ax*By-Ay*Bx)
    poly=sp.Poly(J, x,y)
    for coeff in poly.as_dict().values():
        eqs.append(coeff)
    try:
        sol = sp.solve([e for e in eqs if e!=0], list(a_syms)+list(b_syms), dict=True)
        if sol:
            print(f"Found solution for SA={SA}, SB={SB}: {sol}")
    except Exception as e:
        continue

print("Done sparse enum d3")
