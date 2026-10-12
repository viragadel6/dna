import sympy as sp, random

x,y=sp.symbols('x y')

a20,a11,a02,b20, a30,a21,a12,a03,b30 = sp.symbols('a20 a11 a02 b20 a30 a21 a12 a03 b30')

b11 = -2*a20
b02 = -a11/2
A2 = a20*x**2 + a11*x*y + a02*y**2
B2 = b20*x**2 + b11*x*y + b02*y**2

A3 = a30*x**3 + a21*x**2*y + a12*x*y**2 + a03*y**3

JacA2B2 = sp.diff(A2,x)*sp.diff(B2,y) - sp.diff(A2,y)*sp.diff(B2,x)
polyJac = sp.Poly(JacA2B2, x,y)
c_x2 = polyJac.coeff_monomial(x**2)
c_xy = polyJac.coeff_monomial(x*y)
c_y2 = polyJac.coeff_monomial(y**2)

b21_expr = -3*a30 - c_x2
b12_expr = -a21 - c_xy/2
b03_expr = (-a12 - c_y2)/3

B3 = b30*x**3 + b21_expr*x**2*y + b12_expr*x*y**2 + b03_expr*y**3

JacA2B3 = sp.diff(A2,x)*sp.diff(B3,y) - sp.diff(A2,y)*sp.diff(B3,x)
JacA3B2 = sp.diff(A3,x)*sp.diff(B2,y) - sp.diff(A3,y)*sp.diff(B2,x)
J3_rhs = sp.expand(JacA2B3 + JacA3B2)
polyJ3 = sp.Poly(J3_rhs, x,y)
c_x3 = polyJ3.coeff_monomial(x**3)
c_x2y = polyJ3.coeff_monomial(x**2*y)
c_xy2 = polyJ3.coeff_monomial(x*y**2)
c_y3 = polyJ3.coeff_monomial(y**3)

b31_expr = -c_x3
b22_expr = -c_x2y/2
b13_expr = -c_xy2/3
b04_expr = -c_y3/4

# b40 from collision: B(1,0)=b20+b30+b40=0 => b40 = -b20 -b30
b40_expr = -b20 - b30

B4 = b40_expr*x**4 + b31_expr*x**3*y + b22_expr*x**2*y**2 + b13_expr*x*y**3 + b04_expr*y**4

JacA2B4 = sp.diff(A2,x)*sp.diff(B4,y) - sp.diff(A2,y)*sp.diff(B4,x)
JacA3B3 = sp.diff(A3,x)*sp.diff(B3,y) - sp.diff(A3,y)*sp.diff(B3,x)
J4 = sp.expand(JacA2B4 + JacA3B3)
polyJ4 = sp.Poly(J4, x,y)
eqsJ4 = list(polyJ4.as_dict().values())
print(f"J4 eqs: {eqsJ4}")

JacA3B4 = sp.diff(A3,x)*sp.diff(B4,y) - sp.diff(A3,y)*sp.diff(B4,x)
J5 = sp.expand(JacA3B4)
polyJ5 = sp.Poly(J5, x,y)
eqsJ5 = list(polyJ5.as_dict().values())
print(f"J5 eqs: {eqsJ5}")

A_coll = (A2 + A3).subs({x:1,y:0})
eq_coll1 = A_coll + 1
print(f"coll1: {eq_coll1}")

all_eqs = eqsJ4 + eqsJ5 + [eq_coll1]
print(f"Total eqs: {len(all_eqs)}")
vars_list = [a20,a11,a02,b20, a30,a21,a12,a03,b30]
print(f"Vars: {len(vars_list)}")

# Try nsolve with 9 eqs out of 12
# Pick J4 (5) + J5 first 3 (3) + coll1 (1) =9
eqs_square = eqsJ4 + eqsJ5[:3] + [eq_coll1]
print(f"Square eqs: {len(eqs_square)}")

for attempt in range(20):
    print(f"\nAttempt {attempt}")
    init = [random.uniform(-2,2) for _ in range(len(vars_list))]
    try:
        sol = sp.nsolve(eqs_square, vars_list, init, tol=1e-14, maxsteps=300, prec=80)
        print("Found solution:", sol)
        subs = {v:val for v,val in zip(vars_list, sol)}
        # Check remaining eqs
        rem = eqsJ5[3:]
        residuals = [complex(eq.subs(subs).evalf()) for eq in rem]
        print("Remaining residuals:", residuals)
        if all(abs(r) < 1e-6 for r in residuals):
            print("*** FULL SOLUTION FOUND ***")
            A_full = (A2 + A3).subs(subs)
            B_full = (B2 + B3 + B4).subs(subs)
            print("A:", sp.expand(A_full))
            print("B:", sp.expand(B_full))
            Ax = sp.diff(A_full,x)
            Ay = sp.diff(A_full,y)
            Bx = sp.diff(B_full,x)
            By = sp.diff(B_full,y)
            J = sp.expand(Ax+By+Ax*By-Ay*Bx)
            print("J:", J)
            print("A(1,0):", A_full.subs({x:1,y:0}))
            print("B(1,0):", B_full.subs({x:1,y:0}))
            break
    except Exception as e:
        print("nsolve failed:", e)
        continue

print("\nDone")
