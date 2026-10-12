import sympy as sp

x,y=sp.symbols('x y')

# p0,p1,p2,q0,q1,q2 as polynomials in y degree <=2
# p0 = a02 y^2
# p1 = a11 y + a12 y^2
# p2 = -1 + a21 y + a22 y^2  (a20=-1)
# q0 = b02 y^2
# q1 = b11 y + b12 y^2
# q2 = 0 + b21 y + b22 y^2 (b20=0)

a02,a11,a12,a21,a22, b02,b11,b12,b21,b22 = sp.symbols('a02 a11 a12 a21 a22 b02 b11 b12 b21 b22')

p0 = a02*y**2
p1 = a11*y + a12*y**2
p2 = -1 + a21*y + a22*y**2

q0 = b02*y**2
q1 = b11*y + b12*y**2
q2 = b21*y + b22*y**2

A = p0 + x*p1 + x**2*p2
B = q0 + x*q1 + x**2*q2

print("A:", sp.expand(A))
print("B:", sp.expand(B))
print("A(1,0):", A.subs({x:1,y:0}))
print("B(1,0):", B.subs({x:1,y:0}))

Ax=sp.diff(A,x)
Ay=sp.diff(A,y)
Bx=sp.diff(B,x)
By=sp.diff(B,y)
J=sp.expand(Ax+By+Ax*By-Ay*Bx)
print("\nJ:", J)
poly=sp.Poly(J, x,y)
print("J coeffs:", poly.as_dict())
eqs=list(poly.as_dict().values())
print(f"Number eqs: {len(eqs)}")
# Solve
sols=sp.solve(eqs, [a02,a11,a12,a21,a22,b02,b11,b12,b21,b22], dict=True)
print("Solutions:", sols)

# Try degree up to 3 in y
print("\n\n--- Degree up to 3 in y ---")
a02,a03,a11,a12,a13,a21,a22,a23, b02,b03,b11,b12,b13,b21,b22,b23 = sp.symbols('a02 a03 a11 a12 a13 a21 a22 a23 b02 b03 b11 b12 b13 b21 b22 b23')
p0 = a02*y**2 + a03*y**3
p1 = a11*y + a12*y**2 + a13*y**3
p2 = -1 + a21*y + a22*y**2 + a23*y**3
q0 = b02*y**2 + b03*y**3
q1 = b11*y + b12*y**2 + b13*y**3
q2 = b21*y + b22*y**2 + b23*y**3

A = p0 + x*p1 + x**2*p2
B = q0 + x*q1 + x**2*q2

Ax=sp.diff(A,x)
Ay=sp.diff(A,y)
Bx=sp.diff(B,x)
By=sp.diff(B,y)
J=sp.expand(Ax+By+Ax*By-Ay*Bx)
poly=sp.Poly(J, x,y)
eqs=list(poly.as_dict().values())
print(f"Number eqs: {len(eqs)}, vars: 16")
# Try solve with Groebner? Might be heavy
# Let's try to solve with sympy's solve for subset
# Use nsolve with random initial guesses for numeric search
import random
vars_list = [a02,a03,a11,a12,a13,a21,a22,a23,b02,b03,b11,b12,b13,b21,b22,b23]
# Pick 16 eqs (need square)
# J has degree up to? x^2 * y^3 etc, total monomials maybe many, but we have 16 vars, need 16 eqs
# Let's take first 16 eqs
eqs_square = eqs[:16]
print(f"Square eqs: {len(eqs_square)}")
for attempt in range(10):
    init = [random.uniform(-2,2) for _ in range(len(vars_list))]
    try:
        sol = sp.nsolve(eqs_square, vars_list, init, tol=1e-14, maxsteps=200, prec=50)
        print(f"Attempt {attempt} found solution: {sol}")
        subs = {v:val for v,val in zip(vars_list, sol)}
        # Check remaining eqs
        rem = eqs[16:]
        residuals = [complex(eq.subs(subs).evalf()) for eq in rem]
        print("Remaining residuals:", residuals)
        if all(abs(r)<1e-6 for r in residuals):
            print("*** FULL SOLUTION ***")
            A_full = A.subs(subs)
            B_full = B.subs(subs)
            print("A:", sp.expand(A_full))
            print("B:", sp.expand(B_full))
            break
    except Exception as e:
        print(f"Attempt {attempt} nsolve failed: {e}")
        continue

print("Done")
