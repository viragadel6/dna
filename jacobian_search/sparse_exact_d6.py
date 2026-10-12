import sympy as sp, random, itertools

x,y=sp.symbols('x y')

monoms_d6 = []
for i in range(0,7):
    for j in range(0,7):
        if 2 <= i+j <=6:
            monoms_d6.append((i,j))
print(f"Monoms d6: {len(monoms_d6)}")

def monom_expr(i,j):
    return x**i * y**j

def try_support(SA, SB):
    a_syms = sp.symbols(f'a0:{len(SA)}')
    b_syms = sp.symbols(f'b0:{len(SB)}')
    A = sum(a_syms[k]*monom_expr(*SA[k]) for k in range(len(SA)))
    B = sum(b_syms[k]*monom_expr(*SB[k]) for k in range(len(SB)))
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
    try:
        sol = sp.solve([e for e in eqs if e!=0], list(a_syms)+list(b_syms), dict=True)
        if sol:
            print(f"Found solution SA={SA}, SB={SB}: {sol}")
            return True
    except:
        pass
    return False

random.seed(1)
# Generate all supports size 5: C(27,5)=80751, too many, random sample
for trial in range(50):
    SA = tuple(random.sample(monoms_d6, 5))
    SB = tuple(random.sample(monoms_d6, 5))
    if not any(j==0 for i,j in SA):
        continue
    if not any(j==0 for i,j in SB):
        continue
    print(f"Trial {trial}: SA={SA}, SB={SB}")
    found = try_support(SA, SB)
    if found:
        break

print("Done d6 sparse")
