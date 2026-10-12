import sympy as sp
x,y,lam,t = sp.symbols('x y lam t')
s = lam*t+1
# Avoid division: multiply? Keep rational
M = x - y/lam
L = x + t*y
A2 = -1/s * M**2
B2 = -lam/s * M**2
c = -lam*t/s
cp = lam/s
A3 = c * L**3
B3 = cp * L**3
A = A2 + A3
B = B2 + B3
Ax = sp.diff(A,x)
Ay = sp.diff(A,y)
Bx = sp.diff(B,x)
By = sp.diff(B,y)
J = sp.expand(Ax + By + Ax*By - Ay*Bx)
print("J:", J)
# Simplify: factor s etc
J_simpl = sp.simplify(J)
print("J simplified:", J_simpl)
# Expand numerator
# Put over common denominator? Let's compute with sympy rational: A,B are rational in lam,t but polynomial in x,y with rational coeffs.
# Compute polynomial in x,y with coefficients rational in lam,t
poly = sp.Poly(J, x, y)
print("coeffs:", poly.as_dict())
# Check if identically zero for all lam,t?
# Set lam=1,t=0 as example
lam_val=1
t_val=0
A_test = A.subs({lam:lam_val, t:t_val})
B_test = B.subs({lam:lam_val, t:t_val})
print("\nTest lam=1,t=0:")
print("A:", sp.expand(A_test))
print("B:", sp.expand(B_test))
Ax_test = sp.diff(A_test,x)
Ay_test = sp.diff(A_test,y)
Bx_test = sp.diff(B_test,x)
By_test = sp.diff(B_test,y)
J_test = sp.expand(Ax_test+By_test+Ax_test*By_test-Ay_test*Bx_test)
print("J_test:", J_test)
# Collision
print("A(1,0):", A_test.subs({x:1,y:0}))
print("B(1,0):", B_test.subs({x:1,y:0}))
# F1 = x+A, F2=y+B
F1_test = x + A_test
F2_test = y + B_test
print("F1:", sp.expand(F1_test))
print("F2:", sp.expand(F2_test))
print("F(0,0):", F1_test.subs({x:0,y:0}), F2_test.subs({x:0,y:0}))
print("F(1,0):", F1_test.subs({x:1,y:0}), F2_test.subs({x:1,y:0}))

# Now general J: maybe we missed J3,J4? Let's compute J3,J4 symbolic general lam,t
# Compute Jacs
Jac_A2_B3 = sp.diff(A2,x)*sp.diff(B3,y) - sp.diff(A2,y)*sp.diff(B3,x)
Jac_A3_B2 = sp.diff(A3,x)*sp.diff(B2,y) - sp.diff(A3,y)*sp.diff(B2,x)
J3 = sp.expand(Jac_A2_B3 + Jac_A3_B2)
print("\nJ3 general:", sp.simplify(J3))
Jac_A3_B3 = sp.diff(A3,x)*sp.diff(B3,y) - sp.diff(A3,y)*sp.diff(B3,x)
print("J4 general (should be zero):", sp.simplify(Jac_A3_B3))

# So J = J2 + J3 + J4? Actually J = J1+J2+J3+J4, J1=0, J2=0 by construction, so J = J3+J4. Since J4=0, J=J3. So J3 must be zero for Keller.
# Let's see J3 expression factored
print("\nJ3 factored:", sp.factor(J3))
# Solve for lam,t making J3 zero
# J3 is homogeneous degree3: coefficients must be zero
polyJ3 = sp.Poly(J3, x, y)
print("J3 coeffs:", polyJ3.as_dict())
# Set each coeff zero, solve for lam,t?
# Coeffs are rational functions in lam,t
# Let's get numerator
for monom, coeff in polyJ3.as_dict().items():
    num = sp.simplify(sp.together(coeff).as_numer_denom()[0])
    den = sp.simplify(sp.together(coeff).as_numer_denom()[1])
    print(monom, "coeff num:", sp.factor(num), "den:", den)

# So we need numerators zero. Let's attempt solve for t with generic lam
# Numerators likely have factor?
# Let's attempt sympy solve for t
coeffs = list(polyJ3.as_dict().values())
# Multiply by s^2 etc to clear denominators
# Let's clear denominators: s = lam*t+1, lam denominator
# Multiply J3 by s^2 * lam^? Let's compute J3 * s^2 * lam
J3_scaled = sp.expand(J3 * s**2 * lam)
print("\nJ3_scaled:", J3_scaled)
poly_scaled = sp.Poly(J3_scaled, x,y)
print(poly_scaled.as_dict())

# For lam=1, solve for t making J3 zero
lam_val=1
J3_lam1 = J3.subs({lam:lam_val})
print("\nJ3 lam=1:", J3_lam1)
poly_lam1 = sp.Poly(J3_lam1, x,y)
for monom, coeff in poly_lam1.as_dict().items():
    print(monom, sp.factor(coeff))

# So for lam=1, J3 =? Let's compute with t variable
J3_lam1_t = sp.expand(J3_lam1)
# Solve coeff zero
# Coeffs are functions of t
# Extract equations for t
eqs_t = []
for coeff in poly_lam1.as_dict().values():
    eqs_t.append(sp.simplify(coeff))
print("eqs_t:", eqs_t)
# Solve
sol_t = sp.solve(eqs_t, [t], dict=True)
print("sol_t lam=1:", sol_t)
# Might need solve for t numeric
# Use resultant?
# Let's try numeric: find t such that all coeffs zero
# The system is overdetermined, maybe only t=0 works? Let's test t=0
print("J3 lam=1 t=0:", J3_lam1.subs({t:0}))
# t = -1/lam gives s=0 singular, not allowed

# Try general lam, solve for t
# From J3_scaled coeffs: let's get equations
eqs_general = list(poly_scaled.as_dict().values())
print("\nGeneral eqs for J3_scaled zero:")
for e in eqs_general:
    print(sp.factor(e))

# Solve for t in terms of lam?
# Let's attempt to solve one equation for t
# Take first equation
e0 = eqs_general[0]
print("\ne0:", e0)
# Solve e0=0 for t
sol_e0 = sp.solve(e0, t)
print("sol_e0:", sol_e0)

# Let's brute force search for rational lam,t that make J3 zero: maybe lam arbitrary, t=0 is solution?
# Check t=0 generally
J3_t0 = J3.subs({t:0})
print("\nJ3 t=0 general lam:", sp.simplify(J3_t0))
poly_t0 = sp.Poly(J3_t0, x,y)
print(poly_t0.as_dict())
# t=0 => J3 =?
# For t=0, s=1, M=x - y/lam, L=x
# A2=-M^2, B2=-lam M^2, A3=0, B3=lam x^3
# Then Jac(A2,B3)= Jac(-M^2, lam x^3) = -lam * Jac(M^2, x^3)
# Jac(M^2, x^3)=? Compute: M = x - y/lam, M^2 = x^2 -2xy/lam + y^2/lam^2
# M^2_x=2x-2y/lam, M^2_y=-2x/lam+2y/lam^2
# x^3_x=3x^2, x^3_y=0 => Jac= (2x-2y/lam)*0 - (-2x/lam+2y/lam^2)*3x^2 =3x^2(2x/lam -2y/lam^2)=6 x^3/lam -6 x^2 y /lam^2
# So Jac(A2,B3)= -lam * (6 x^3/lam -6 x^2 y/lam^2)= -6 x^3 +6 x^2 y /lam
# Jac(A3,B2)=0 since A3=0 => J3 = -6 x^3 +6 x^2 y/lam ≠0 generally! So t=0 not solution.

# So our family fails J3.

# Conclusion: d=3 full Keller with collision may have no solution, but we found intermediate family satisfying J1,J2 only.

# Need to search more generally for d=3 without assuming proportionality of top? But we proved top must be proportional, and then J3 imposes extra.

# Let's do exhaustive Groebner for d=3 remaining 11 equations (9 unknowns). Use sympy groebner with elimination.

print("\n--- Attempt full Groebner for d=3 remaining system ---")
# Variables: a20,a11,a02,b20,a30,a21,a12,a03,b30 (9 vars)
a20_,a11_,a02_,b20_,a30_,a21_,a12_,a03_,b30_ = sp.symbols('a20_ a11_ a02_ b20_ a30_ a21_ a12_ a03_ b30_')
# Define expressions for b21,b12,b03 as before
c_x2 = -2*a11_*b20_ -4*a20_**2
c_xy = -4*a02_*b20_ -2*a11_*a20_
c_y2 = 4*a02_*a20_ - a11_**2
b21_ = -3*a30_ - c_x2
b12_ = -a21_ - c_xy/2
b03_ = (-a12_ - c_y2)/3

A2_ = a20_*x**2 + a11_*x*y + a02_*y**2
B2_ = b20_*x**2 -2*a20_*x*y - a11_/2*y**2
A3_ = a30_*x**3 + a21_*x**2*y + a12_*x*y**2 + a03_*y**3
B3_ = b30_*x**3 + b21_*x**2*y + b12_*x*y**2 + b03_*y**3

Jac_A2_B3_ = sp.diff(A2_,x)*sp.diff(B3_,y) - sp.diff(A2_,y)*sp.diff(B3_,x)
Jac_A3_B2_ = sp.diff(A3_,x)*sp.diff(B2_,y) - sp.diff(A3_,y)*sp.diff(B2_,x)
J3_ = sp.expand(Jac_A2_B3_ + Jac_A3_B2_)
Jac_A3_B3_ = sp.diff(A3_,x)*sp.diff(B3_,y) - sp.diff(A3_,y)*sp.diff(B3_,x)
J4_ = sp.expand(Jac_A3_B3_)

polyJ3_ = sp.Poly(J3_, x,y)
eqsJ3_ = list(polyJ3_.as_dict().values())
polyJ4_ = sp.Poly(J4_, x,y)
eqsJ4_ = list(polyJ4_.as_dict().values())

eq_coll1_ = a20_ + a30_ +1
eq_coll2_ = b20_ + b30_

all_eqs_ = eqsJ3_ + eqsJ4_ + [eq_coll1_, eq_coll2_]
print(f"Total eqs: {len(all_eqs_)}")
for e in all_eqs_:
    print(e)

# Try Groebner basis
try:
    G = sp.groebner(all_eqs_, a20_,a11_,a02_,b20_,a30_,a21_,a12_,a03_,b30_, order='lex')
    print("Groebner basis:", G)
    if 1 in G.polys:
        print("System inconsistent (Groebner contains 1) => no solution for d=3")
    else:
        print("Groebner computed, not containing 1, may have solutions")
except Exception as e:
    print("Groebner failed:", e)

# Try solving with elimination via resultant or using sympy's solve with dict?
# Could try to solve J4 first: Jac(A3,B3)=0 => proportional condition already analyzed, but let's keep.

print("\nDone")
