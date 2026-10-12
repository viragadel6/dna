import sympy as sp
import random
import cmath

x,y = sp.symbols('x y')

# Degree 4: coefficients
# A2: a20,a11,a02
# B2: b20,b11,b02 with constraints b11=-2a20, b02=-a11/2 from J1
# Let's keep free: a20,a11,a02,b20
# A3: a30,a21,a12,a03
# B3: b30,b21,b12,b03 with b21,b12,b03 expressed via J2
# A4: a40,a31,a22,a13,a04
# B4: b40,b31,b22,b13,b04 with linear relations from J3

# Let's derive symbolic expressions for b21,b12,b03 and for B4 linear relations

# Step1: A2,B2
a20,a11,a02,b20 = sp.symbols('a20 a11 a02 b20')
b11 = -2*a20
b02 = -a11/2

A2 = a20*x**2 + a11*x*y + a02*y**2
B2 = b20*x**2 + b11*x*y + b02*y**2

# A3,B3 generic
a30,a21,a12,a03,b30 = sp.symbols('a30 a21 a12 a03 b30')
a21_,a12_,a03_,b30_ = a21,a12,a03,b30
# Actually we need symbols for A3 coeffs
# A3 = a30 x^3 + a21 x^2 y + a12 x y^2 + a03 y^3
# B3 = b30 x^3 + b21 x^2 y + b12 x y^2 + b03 y^3
b21,b12,b03 = sp.symbols('b21 b12 b03')
A3 = a30*x**3 + a21*x**2*y + a12*x*y**2 + a03*y**3
B3 = b30*x**3 + b21*x**2*y + b12*x*y**2 + b03*y**3

# J2: (A3)_x + (B3)_y + Jac(A2,B2)=0
JacA2B2 = sp.diff(A2,x)*sp.diff(B2,y) - sp.diff(A2,y)*sp.diff(B2,x)
J2 = sp.expand(sp.diff(A3,x) + sp.diff(B3,y) + JacA2B2)
polyJ2 = sp.Poly(J2, x,y)
# coeffs for x^2, xy, y^2
# Solve for b21,b12,b03
# J2 = (3a30 + b21 + coeff_x2_Jac) x^2 + (2a21+2b12+coeff_xy) xy + (a12+3b03+coeff_y2) y^2
coeffs_Jac = sp.Poly(JacA2B2, x,y).as_dict()
c_x2 = coeffs_Jac.get((2,0),0)
c_xy = coeffs_Jac.get((1,1),0)
c_y2 = coeffs_Jac.get((0,2),0)
print("JacA2B2 coeffs:", c_x2,c_xy,c_y2)

# So b21 = -3a30 - c_x2, etc.
b21_expr = -3*a30 - c_x2
b12_expr = -a21 - c_xy/2
b03_expr = (-a12 - c_y2)/3

print("b21,b12,b03 expr:", b21_expr, b12_expr, b03_expr)

B3_sub = B3.subs({b21:b21_expr, b12:b12_expr, b03:b03_expr})

# Now J3: (A4)_x + (B4)_y + Jac(A2,B3)+Jac(A3,B2)=0
# A4 generic
a40,a31,a22,a13,a04 = sp.symbols('a40 a31 a22 a13 a04')
b40,b31,b22,b13,b04 = sp.symbols('b40 b31 b22 b13 b04')
A4 = a40*x**4 + a31*x**3*y + a22*x**2*y**2 + a13*x*y**3 + a04*y**4
B4 = b40*x**4 + b31*x**3*y + b22*x**2*y**2 + b13*x*y**3 + b04*y**4

JacA2B3 = sp.diff(A2,x)*sp.diff(B3_sub,y) - sp.diff(A2,y)*sp.diff(B3_sub,x)
JacA3B2 = sp.diff(A3,x)*sp.diff(B2,y) - sp.diff(A3,y)*sp.diff(B2,x)
J3_rhs = sp.expand(JacA2B3 + JacA3B2)
# J3 equation: (A4)_x + (B4)_y = -J3_rhs
A4x = sp.diff(A4,x)
B4y = sp.diff(B4,y)
J3 = sp.expand(A4x + B4y + J3_rhs)
polyJ3 = sp.Poly(J3, x,y)
print("\nJ3 poly dict:", polyJ3.as_dict())
# J3 is degree3 homogeneous: monomials x^3,x^2y,xy^2,y^3 (4 eqs)
# Solve for B4 coefficients in terms of A4 and known
# A4x = 4a40 x^3 +3a31 x^2 y +2a22 x y^2 + a13 y^3
# B4y = b31 x^3 +2b22 x^2 y +3b13 x y^2 +4b04 y^3
# So:
# x^3: 4a40 + b31 + coeff_J3_x3 =0 => b31 = -4a40 - coeff_x3
# x^2y: 3a31 +2b22 + coeff_x2y=0 => b22 = (-3a31 - coeff_x2y)/2
# xy^2: 2a22 +3b13 + coeff_xy2=0 => b13 = (-2a22 - coeff_xy2)/3
# y^3: a13 +4b04 + coeff_y3=0 => b04 = (-a13 - coeff_y3)/4

# Get coeffs of J3_rhs
coeffs_J3_rhs = sp.Poly(J3_rhs, x,y).as_dict()
c_x3 = coeffs_J3_rhs.get((3,0),0)
c_x2y = coeffs_J3_rhs.get((2,1),0)
c_xy2 = coeffs_J3_rhs.get((1,2),0)
c_y3 = coeffs_J3_rhs.get((0,3),0)
print("J3_rhs coeffs:", c_x3,c_x2y,c_xy2,c_y3)

b31_expr = -4*a40 - c_x3
b22_expr = (-3*a31 - c_x2y)/2
b13_expr = (-2*a22 - c_xy2)/3
b04_expr = (-a13 - c_y3)/4

B4_sub = B4.subs({b31:b31_expr, b22:b22_expr, b13:b13_expr, b04:b04_expr})

# Now remaining equations: J4,J5,J6
# J4: Jac(A2,B4)+Jac(A3,B3)+Jac(A4,B2)=0 degree4 (5 eqs)
# J5: Jac(A3,B4)+Jac(A4,B3)=0 degree5 (6 eqs)
# J6: Jac(A4,B4)=0 degree6 (7 eqs)
# Total 18 eqs

JacA2B4 = sp.diff(A2,x)*sp.diff(B4_sub,y) - sp.diff(A2,y)*sp.diff(B4_sub,x)
JacA3B3 = sp.diff(A3,x)*sp.diff(B3_sub,y) - sp.diff(A3,y)*sp.diff(B3_sub,x)
JacA4B2 = sp.diff(A4,x)*sp.diff(B2,y) - sp.diff(A4,y)*sp.diff(B2,x)
J4 = sp.expand(JacA2B4 + JacA3B3 + JacA4B2)
polyJ4 = sp.Poly(J4, x,y)
eqsJ4 = list(polyJ4.as_dict().values())
print(f"\nJ4 has {len(eqsJ4)} eqs")

JacA3B4 = sp.diff(A3,x)*sp.diff(B4_sub,y) - sp.diff(A3,y)*sp.diff(B4_sub,x)
JacA4B3 = sp.diff(A4,x)*sp.diff(B3_sub,y) - sp.diff(A4,y)*sp.diff(B3_sub,x)
J5 = sp.expand(JacA3B4 + JacA4B3)
polyJ5 = sp.Poly(J5, x,y)
eqsJ5 = list(polyJ5.as_dict().values())
print(f"J5 has {len(eqsJ5)} eqs")

JacA4B4 = sp.diff(A4,x)*sp.diff(B4_sub,y) - sp.diff(A4,y)*sp.diff(B4_sub,x)
J6 = sp.expand(JacA4B4)
polyJ6 = sp.Poly(J6, x,y)
eqsJ6 = list(polyJ6.as_dict().values())
print(f"J6 has {len(eqsJ6)} eqs")

# Collision
A_coll = (A2 + A3 + A4).subs({x:1,y:0})  # = a20 + a30 + a40
B_coll = (B2 + B3_sub + B4_sub).subs({x:1,y:0})  # = b20 + b30 + b40
print(f"\nA_coll: {A_coll}, B_coll: {B_coll}")
eq_coll1 = A_coll + 1
eq_coll2 = B_coll

all_eqs = eqsJ4 + eqsJ5 + eqsJ6 + [eq_coll1, eq_coll2]
print(f"\nTotal remaining eqs: {len(all_eqs)}")
# Variables: a20,a11,a02,b20, a30,a21,a12,a03,b30, a40,a31,a22,a13,a04,b40
# That's 15 variables
vars_list = [a20,a11,a02,b20, a30,a21,a12,a03,b30, a40,a31,a22,a13,a04,b40]
print(f"Variables: {len(vars_list)}")

# Now attempt numeric solving with random initial guesses using nsolve
# For nsolve we need square system: choose 15 equations out of 20
# Let's pick: J4 (5) + J5 (6) + J6 first 2? + collision 2 =15
# But J6 has 7 eqs, maybe pick 4 of them? Let's pick J4 (5) + J5 (6) + first 2 of J6 (2) + coll (2) =15

eqs_for_numerics = eqsJ4 + eqsJ5[:6] + eqsJ6[:2] + [eq_coll1, eq_coll2]
print(f"eqs_for_numerics len: {len(eqs_for_numerics)}")

# Try nsolve with random initial complex guesses
# Use sympy nsolve

# Prepare for nsolve: need to convert to numerical
# Let's attempt with real initial guesses

for attempt in range(5):
    print(f"\n--- Attempt {attempt} nsolve ---")
    # random initial guess
    init_guess = [random.uniform(-2,2) for _ in range(len(vars_list))]
    # Use complex? nsolve works with real by default, but we can allow complex by using complex initial?
    # Try real first
    try:
        sol = sp.nsolve(eqs_for_numerics, vars_list, init_guess, tol=1e-14, maxsteps=200, prec=50)
        print("Found solution (real):", sol)
        # Check remaining eqs
        subs_dict = {var: val for var, val in zip(vars_list, sol)}
        # Evaluate remaining J6 eqs
        remaining = eqsJ6[2:]
        res_remaining = [abs(complex(eq.subs(subs_dict).evalf())) for eq in remaining]
        print("Remaining J6 residuals:", res_remaining)
        # Also check full J?
        # If residuals small, we have full solution!
        if all(r < 1e-6 for r in res_remaining):
            print("*** FOUND FULL SOLUTION ***")
            # Compute A,B
            A_full = (A2 + A3 + A4).subs(subs_dict)
            B_full = (B2 + B3_sub + B4_sub).subs(subs_dict)
            print("A:", sp.expand(A_full))
            print("B:", sp.expand(B_full))
            # Check Jacobian
            Ax = sp.diff(A_full,x)
            Ay = sp.diff(A_full,y)
            Bx = sp.diff(B_full,x)
            By = sp.diff(B_full,y)
            J_full = sp.expand(Ax+By+Ax*By-Ay*Bx)
            print("J residual:", J_full)
            # Check collision
            print("A(1,0):", A_full.subs({x:1,y:0}))
            print("B(1,0):", B_full.subs({x:1,y:0}))
            break
    except Exception as e:
        print(f"nsolve failed: {e}")
        continue

print("\nDone numeric attempts")
