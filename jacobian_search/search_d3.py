import sympy as sp

# Variables
x,y = sp.symbols('x y')

# Coefficients for degree 2 and 3
# A2 = a20 x^2 + a11 xy + a02 y^2
# B2 = b20 x^2 + b11 xy + b02 y^2
a20,a11,a02,b20,b11,b02 = sp.symbols('a20 a11 a02 b20 b11 b02')
# A3 = a30 x^3 + a21 x^2 y + a12 x y^2 + a03 y^3
# B3 = b30 x^3 + b21 x^2 y + b12 x y^2 + b03 y^3
a30,a21,a12,a03,b30,b21,b12,b03 = sp.symbols('a30 a21 a12 a03 b30 b21 b12 b03')

A2 = a20*x**2 + a11*x*y + a02*y**2
B2 = b20*x**2 + b11*x*y + b02*y**2

A3 = a30*x**3 + a21*x**2*y + a12*x*y**2 + a03*y**3
B3 = b30*x**3 + b21*x**2*y + b12*x*y**2 + b03*y**3

A = A2 + A3
B = B2 + B3

# Jacobian condition: A_x + B_y + A_x B_y - A_y B_x =0
Ax = sp.diff(A, x)
Ay = sp.diff(A, y)
Bx = sp.diff(B, x)
By = sp.diff(B, y)

J = Ax + By + Ax*By - Ay*Bx
J = sp.expand(J)
poly = sp.Poly(J, x, y)
print("J poly terms:", poly.terms())
coeffs = poly.coeffs()
# Get equations: all coefficients zero
eqs = []
for monom, coeff in poly.as_dict().items():
    eqs.append(coeff)

print(f"Number equations from J: {len(eqs)}")
for i, e in enumerate(eqs):
    print(f"eq{i}: {e}")

# Collision
eq_coll1 = (A.subs({x:1,y:0}) + 1)  # =0
eq_coll2 = B.subs({x:1,y:0})  # =0
print(f"coll1: {eq_coll1}")
print(f"coll2: {eq_coll2}")
eqs.append(eq_coll1)
eqs.append(eq_coll2)

print("\nTotal eqs:", len(eqs))
print("Variables:", 14)

# Try to solve stepwise using homogeneous decomposition
# J1: degree 1
# Extract degree 1 part: Ax from A2 is linear, By from B2 linear
J1 = sp.expand(sp.diff(A2,x) + sp.diff(B2,y))
print("\nJ1:", J1)
# J1 = (2 a20 + b11) x + (a11 + 2 b02) y
# So equations:
# 2 a20 + b11 =0
# a11 + 2 b02 =0

# Solve these for b11, b02
sol_b11 = -2*a20
sol_b02 = -a11/2
print(f"b11 = {sol_b11}, b02 = {sol_b02}")

# Substitute
subs1 = {b11: sol_b11, b02: sol_b02}
A2s = A2
B2s = B2.subs(subs1)
print("B2s:", B2s)

# J2: degree 2 part
# (A3)_x + (B3)_y + Jac(A2,B2)
A3x = sp.diff(A3,x)
B3y = sp.diff(B3,y)
JacA2B2 = sp.diff(A2s,x)*sp.diff(B2s,y) - sp.diff(A2s,y)*sp.diff(B2s,x)
J2 = sp.expand(A3x + B3y + JacA2B2)
print("\nJ2:", J2)
polyJ2 = sp.Poly(J2, x, y)
eqsJ2 = list(polyJ2.as_dict().values())
print("eqsJ2:", eqsJ2)

# These are linear in A3,B3 coeffs: 3 equations, 8 unknowns
# Solve linear system
# J2 = coeff x^2 + coeff xy + coeff y^2
# Extract
# Write as linear equations
# Use sympy linear_eq_to_matrix
vars_J2 = [a30,a21,a12,a03,b30,b21,b12,b03]
M, b = sp.linear_eq_to_matrix(eqsJ2, vars_J2)
print("M shape", M.shape)
print(M)
print(b)
# Solve: M * vars = - (remaining part from Jac?) Actually b already includes Jac part? Wait eqsJ2 includes Jac, so equation is expression =0.
# M*vars + const =0 where const from Jac
# So we have M*vars = -const_part
# Let's compute
# Actually linear_eq_to_matrix returns M, b such that M*vars + b =0? Need check
# According docs: returns M,b where M*vars = b ? Let's test: linear_eq_to_matrix([2*x+y-2], [x,y]) -> M=[[2,1]], b=[2]
# So equation M*vars = b? Actually 2x+y =2 => b=2
# So our eqs are sum =0, so we need M*vars + const_Jac =0 => M*vars = -const_Jac
# Let's compute directly: split Jac part

Jac_poly = sp.Poly(JacA2B2, x, y)
print("JacA2B2:", JacA2B2)
# Its coefficients are quadratic in a20,a11,a02,b20
# The A3x+B3y part: A3x = 3 a30 x^2 +2 a21 xy + a12 y^2, B3y = b21 x^2 +2 b12 xy +3 b03 y^2
# So:
# x^2: 3 a30 + b21 + coeff_Jac_x2 =0
# xy: 2 a21 +2 b12 + coeff_Jac_xy =0
# y^2: a12 +3 b03 + coeff_Jac_y2 =0

# Let's compute Jac coefficients
c_x2 = Jac_poly.coeff_monomial(x**2)
c_xy = Jac_poly.coeff_monomial(x*y)
c_y2 = Jac_poly.coeff_monomial(y**2)
print(f"Jac coeffs: x2={c_x2}, xy={c_xy}, y2={c_y2}")

# So linear solution:
# b21 = -3 a30 - c_x2
# b12 = - a21 - c_xy/2
# b03 = (-a12 - c_y2)/3

# Let's param: keep a30,a21,a12,a03,b30 as free, express b21,b12,b03

b21_expr = -3*a30 - c_x2
b12_expr = -a21 - c_xy/2
b03_expr = (-a12 - c_y2)/3

print(f"b21={b21_expr}\nb12={b12_expr}\nb03={b03_expr}")

# Now remaining equations: J3 degree3 and J4 degree4
# J3 = Jac(A2,B3)+Jac(A3,B2) + (A4)_x+(B4)_y but A4=0 for d=3, so J3 = Jac(A2,B3)+Jac(A3,B2) =0
# J4 = Jac(A3,B3)=0

# Substitute b21,b12,b03 expressions into B3
B3_sub = B3.subs({b21:b21_expr, b12:b12_expr, b03:b03_expr})
# Also B2s already
# Compute Jacs
Jac_A2_B3 = sp.diff(A2s,x)*sp.diff(B3_sub,y) - sp.diff(A2s,y)*sp.diff(B3_sub,x)
Jac_A3_B2 = sp.diff(A3,x)*sp.diff(B2s,y) - sp.diff(A3,y)*sp.diff(B2s,x)
J3 = sp.expand(Jac_A2_B3 + Jac_A3_B2)
print("\nJ3:", J3)
polyJ3 = sp.Poly(J3, x, y)
eqsJ3 = list(polyJ3.as_dict().values())
print(f"J3 has {len(eqsJ3)} equations degree 3 homogeneous (should be 4 monomials x^3,x^2y,xy^2,y^3)")
for e in eqsJ3:
    print(sp.expand(e))

Jac_A3_B3 = sp.diff(A3,x)*sp.diff(B3_sub,y) - sp.diff(A3,y)*sp.diff(B3_sub,x)
J4 = sp.expand(Jac_A3_B3)
print("\nJ4:", J4)
polyJ4 = sp.Poly(J4, x, y)
eqsJ4 = list(polyJ4.as_dict().values())
print(f"J4 has {len(eqsJ4)} eqs")
for e in eqsJ4:
    print(e)

# Now total system after reduction: unknowns = a20,a11,a02,b20, a30,a21,a12,a03,b30
# That's 9 unknowns
# Equations: J3 (4 eqs) + J4 (5 eqs) + collision (2 eqs) =11 eqs
# Overdetermined but nonlinear (degree up to 3?)

# Let's write collision eqs with substitutions
A_coll = (A2s + A3).subs({x:1,y:0})  # = a20 + a30
B_coll = (B2s + B3_sub).subs({x:1,y:0}) # = b20 + b30
print(f"\nCollision: A(1,0) = {A_coll} = -1 => {A_coll}+1=0")
print(f"B(1,0) = {B_coll}=0")

eq_coll1 = A_coll + 1
eq_coll2 = B_coll

all_eqs = eqsJ3 + eqsJ4 + [eq_coll1, eq_coll2]
print(f"\nTotal remaining eqs: {len(all_eqs)}")
# Try Groebner or solving
# Use sympy solve? Might be heavy.

# Let's attempt to compute Groebner basis with lex order for small subset?
# First try to solve with rational elimination using sympy's groebner for part of system?

# For d=3, we suspect no solution. Let's try numeric solving with random?
# Use sympy's solve for some variables?

# Attempt to compute Groebner for J4 alone: Jac(A3,B3)=0
# This implies A3 and B3 are proportional? Let's analyze.
# If Jac(A3,B3)=0, then A3,B3 share common factor structure.
# Since homogeneous degree 3, proportionality condition means there exists lambda such that B3 = lambda A3 if A3 non-zero and irreducible? Actually condition for homogeneous same degree is that they are linearly dependent over C as polynomials? Let's test: Jac(f,g)=0 for homogeneous same degree implies f and g are linearly dependent? Is that true?
# For homogeneous polynomials in 2 variables, f = sum a_i x^{3-i} y^i, g = sum b_i x^{3-i} y^i.
# Jac = f_x g_y - f_y g_x =0.
# This is condition that f,g are algebraically dependent. For binary forms, this implies they have a common factor of degree >=1 and are powers? Let's test: f = x^3, g = x^3 => Jac=0, proportional.
# f = x^2 y, g = x y^2 => Jac = (2xy)(x*? wait compute) => (2xy)(2xy?) Actually f_x=2xy, f_y=x^2, g_x=y^2, g_y=2xy => Jac=4x^2 y^2 - x^2 y^2 =3 x^2 y^2 !=0.
# f = x^3 + y^3, g = (x^3 + y^3) => proportional.
# f = (x^2+xy)^? No.
# Conjecture: for same degree binary forms, Jac=0 iff proportional. Let's test with sympy general condition.

# Compute condition for proportionality: rank of coefficient matrix? Actually if f and g proportional, then there exists lambda: b_i = lambda a_i for all i.
# Does Jac=0 imply this? Let's check dimension: space of pairs (f,g) with Jac=0 has dimension? Proportional pairs have dimension: dim f (4) + 1 (lambda) =5. Total space of pairs is 8 dimensions. Jac=0 gives 5 equations (coeffs of degree 4 homogeneous: 5 monomials). So expected dimension of solution variety is 8-5=3, but proportional variety dimension 5? Something off. Let's compute: Jac degree 4 has 5 coefficients, but they are not independent? There may be dependencies.
# Let's analyze more directly via sympy elimination for general A3,B3.

# For generic A3,B3, Jac condition is 5 equations quadratic.
# Let's try to solve Jac=0 symbolically for A3,B3.

print("\n--- Analyzing Jac(A3,B3)=0 general ---")
# Variables a30,a21,a12,a03,b30,b21,b12,b03 (but we have relations for b21,b12,b03)
# Actually for general A3,B3, without previous substitution, let's analyze
a30_,a21_,a12_,a03_,b30_,b21_,b12_,b03_ = sp.symbols('a30_ a21_ a12_ a03_ b30_ b21_ b12_ b03_')
A3_ = a30_*x**3 + a21_*x**2*y + a12_*x*y**2 + a03_*y**3
B3_ = b30_*x**3 + b21_*x**2*y + b12_*x*y**2 + b03_*y**3
Jac_ = sp.expand(sp.diff(A3_,x)*sp.diff(B3_,y) - sp.diff(A3_,y)*sp.diff(B3_,x))
polyJac_ = sp.Poly(Jac_, x,y)
print(polyJac_.as_dict())

# The system of 5 equations:
# Let's try to find general solution via Groebner with elimination
# Use sympy groebner
# This might be heavy, try with lex order
try:
    G = sp.groebner(list(polyJac_.as_dict().values()), a30_,a21_,a12_,a03_,b30_,b21_,b12_,b03_, order='lex')
    print("Groebner basis for Jac(A3,B3):", G)
except Exception as e:
    print("Groebner failed:", e)

# Alternative: if we set b30 = lambda a30 etc, does it satisfy? Let's test proportionality ansatz: B3 = k * A3
k = sp.symbols('k')
B3_prop = k*A3_
Jac_prop = sp.expand(sp.diff(A3_,x)*sp.diff(B3_prop,y) - sp.diff(A3_,y)*sp.diff(B3_prop,x))
print("Jac for proportional:", sp.simplify(Jac_prop))  # Should be 0

# So proportional is sufficient, but is it necessary?
# Let's search for non-proportional example via solving small case numerically.
# Use random numeric solving with sympy nsolve? Could try.

# Let's attempt to find all solutions where A3 = x^3 (simple) then what B3 must be?
# If A3 = x^3, then A3_x=3x^2, A3_y=0, Jac =3x^2 * B3_y -0 =3x^2 B3_y =0 => B3_y=0 => B3 = b30 x^3 only, proportional.
# If A3 = x^2 y, A3_x=2xy, A3_y=x^2, Jac =2xy B3_y - x^2 B3_x =0 => 2y B3_y - x B3_x =0 => B3 must be homogeneous degree3 satisfying Euler-like. Solve: Let B3 = sum b_i x^{3-i} y^i. Compute B3_x = 3b30 x^2+2b21 xy+ b12 y^2, B3_y = b21 x^2+2 b12 xy+3 b03 y^2. Then 2y B3_y - x B3_x = 2y(b21 x^2+2b12 xy+3b03 y^2) - x(3b30 x^2+2b21 xy+ b12 y^2)= -3b30 x^3 + (2b21-2b21) x^2 y + (4b12 - b12) x y^2 +6 b03 y^3 = -3b30 x^3+3 b12 x y^2+6 b03 y^3=0 => b30=0,b12=0,b03=0, B3 = b21 x^2 y proportional.
# So seems proportional is necessary for these monomials.
# Hypothesis: For binary homogeneous same degree, Jac=0 => proportional. Let's attempt proof: If f,g homogeneous degree d, then Euler: x f_x + y f_y = d f, similarly for g. If Jac=0, then (f_x,g_x) and (f_y,g_y) linearly dependent? Actually Jac=0 means gradient vectors linearly dependent at generic point? That is f_x g_y = f_y g_x. So (f_x,f_y) proportional to (g_x,g_y) as vectors (where defined). So there exists h(x,y) rational such that g_x = h f_x, g_y = h f_y. Then differentiate: etc. Using Euler, d g = x g_x + y g_y = h (x f_x + y f_y)= h d f => g = h f. So h = g/f rational. Then g_x = h f_x + h_x f, but also = h f_x => h_x f=0 => h_x=0 similarly h_y=0 => h constant. So indeed proportional! Provided f not zero divisor. So proof sketch holds for homogeneous same degree. So Jac=0 => proportional.

# Therefore for d=3, B3 = lambda A3 for some lambda (if A3 non-zero). Good.

# Then J3 = Jac(A2,B3)+Jac(A3,B2)= Jac(A2, lambda A3)+Jac(A3,B2)= lambda Jac(A2,A3)+Jac(A3,B2)= Jac(A3, B2 - lambda A2) ??? Wait Jac(A2,lambda A3)= - lambda Jac(A3,A2). So J3 = -lambda Jac(A3,A2)+Jac(A3,B2)= Jac(A3, B2 + lambda? Let's compute: Jac(A3,B2) + lambda Jac(A2,A3) = Jac(A3,B2) - lambda Jac(A3,A2)= Jac(A3, B2 - lambda A2?) No Jac(A3, B2 - lambda A2) = Jac(A3,B2)-lambda Jac(A3,A2). Yes matches. So J3 = Jac(A3, B2 - lambda A2)=0.

# So if A3 non-zero, J3=0 implies B2 - lambda A2 is proportional to A3? But degrees differ: B2 - lambda A2 degree2, A3 degree3. Jac of degree3 and degree2 zero means they are dependent: both powers of same linear form? Let's analyze: If f degree3, g degree2, Jac=0 implies they are algebraically dependent: f = c * h^{3/e}, g = d* h^{2/e} where h homogeneous. So e must divide both 2 and3 => e=1. So h linear: h = u x + v y. Then g = d h^2, f = c h^3. So A3 must be (linear)^3, and B2 - lambda A2 = d' (linear)^2.

# This gives strong structural restriction!

# Let's follow this lemma.

print("\n--- Lemma: Jac homogeneous degrees m,n zero => both powers of same linear form if degrees coprime? ---")

# So for d=3 case:
# If A3 = c L^3, B3 = lambda c L^3 = c' L^3, and B2 - lambda A2 = d L^2 for some linear L.

# Let's attempt to param with L = p x + q y.
# Then A3 = c L^3, B3 = c' L^3
# A2,B2 must satisfy earlier linear constraints and also B2 = lambda A2 + d L^2

# This reduces search drastically!

# Let's implement search with L = x + t y (we can scale) or general p,q.

# Approach: Let L = u x + v y
u,v,c,cp,d,lam = sp.symbols('u v c cp d lam')
# But we also have A2 = a20 x^2 + a11 xy + a02 y^2
# B2 = b20 x^2 -2 a20 xy - a11/2 y^2
# Condition B2 = lam*A2 + d*L^2

# So we have equations for coefficients:
L2 = (u*x+v*y)**2
# Expand and equate coefficients
A2_expr = a20*x**2 + a11*x*y + a02*y**2
B2_expr = b20*x**2 + (-2*a20)*x*y + (-a11/2)*y**2
# B2 - lam*A2 - d*L^2 =0

# So compute polynomial
poly_B2_cond = sp.expand(B2_expr - lam*A2_expr - d*L2)
# Collect coefficients
coeffs_B2_cond = sp.Poly(poly_B2_cond, x,y).as_dict()
print("B2 cond eqs:", coeffs_B2_cond)

# Also A3 = c L^3
L3 = (u*x+v*y)**3
A3_L = c*L3
# Expand to get a30 etc
# A3 = a30 x^3+... So we can equate coefficients to get a30 = c u^3, a21=3c u^2 v, a12=3c u v^2, a03=c v^3
# Similarly B3 = cp L^3

# Now we also have earlier expression for b21,b12,b03 in terms of A2,B2,A3
# But if we are in this structured family, those should be automatically satisfied? Let's check.

# Actually our earlier derivation for b21 etc came from J2 linear equations, which must hold. Let's see if with this ansatz they hold.

# Let's compute J2 condition directly with this ansatz symbolic.

# Let's attempt numeric solving for this restricted family plus collision.

# Variables: a20,a11,a02,b20,u,v,c,cp,d,lam
# Equations: B2 cond (3 eqs) + collision (2 eqs) + maybe J2 automatically? Wait J2 already used to express b21 etc, but if we enforce A3 =c L^3, B3=cp L^3, then J2 condition becomes:
# (A3)_x+(B3)_y + Jac(A2,B2)=0
# Let's compute that condition for general A2,B2,L.

print("\n--- Computing J2 for L-family ---")
# Use generic symbols
a20_s,a11_s,a02_s,b20_s,u_s,v_s,c_s,cp_s,d_s,lam_s = sp.symbols('a20_s a11_s a02_s b20_s u_s v_s c_s cp_s d_s lam_s')
A2_s = a20_s*x**2 + a11_s*x*y + a02_s*y**2
B2_s = b20_s*x**2 -2*a20_s*x*y - a11_s/2*y**2
L_s = u_s*x+v_s*y
A3_s = c_s*L_s**3
B3_s = cp_s*L_s**3
JacA2B2_s = sp.diff(A2_s,x)*sp.diff(B2_s,y) - sp.diff(A2_s,y)*sp.diff(B2_s,x)
J2_s = sp.expand(sp.diff(A3_s,x)+sp.diff(B3_s,y)+JacA2B2_s)
print("J2_s:", J2_s)
polyJ2_s = sp.Poly(J2_s, x,y)
print(polyJ2_s.as_dict())

# This gives 3 equations for unknowns
# Let's attempt to solve these plus B2 cond + collision for a specific choice of u,v? Might be overdetermined but maybe solvable.

# For simplicity, set L = x (v=0) to reduce: L = x
# Then A3 = c x^3, B3 = cp x^3
# Then J2: (A3)_x =3c x^2, (B3)_y=0, Jac(A2,B2) is degree2 homogeneous
# So J2 =3c x^2 + Jac =0 => Jac must be = -3c x^2, i.e., its xy and y^2 coefficients zero, x^2 coeff = -3c

# Let's compute Jac(A2,B2) for general A2,B2:
# A2 = a20 x^2 + a11 xy + a02 y^2
# B2 = b20 x^2 -2a20 xy - a11/2 y^2
# Compute Jac:
# A2_x=2a20 x+ a11 y
# A2_y= a11 x+2a02 y
# B2_x=2b20 x -2a20 y
# B2_y= -2a20 x - a11 y
# Jac = (2a20 x+ a11 y)(-2a20 x - a11 y) - (a11 x+2a02 y)(2b20 x -2a20 y)
# = -(2a20 x+ a11 y)^2 - (a11 x+2a02 y)(2b20 x -2a20 y)
# Expand: 
# = -(4a20^2 x^2+4a20 a11 xy+ a11^2 y^2) - [2a11 b20 x^2 -2a11 a20 xy +4 a02 b20 xy -4 a02 a20 y^2]
# = x^2: -4a20^2 -2a11 b20
# xy: -4a20 a11 +2a11 a20 -4 a02 b20 = -2a20 a11 -4 a02 b20
# y^2: -a11^2 +4 a02 a20
# So Jac = (-4a20^2 -2 a11 b20) x^2 + (-2a20 a11 -4 a02 b20) xy + (-a11^2+4 a02 a20) y^2

# For L=x case, we need xy coeff=0 and y^2 coeff=0, x^2 = -3c

# So system:
# -a11^2+4 a02 a20 =0
# -2a20 a11 -4 a02 b20 =0
# -4a20^2 -2 a11 b20 +3c =0? Actually Jac +3c x^2=0 => -4a20^2 -2a11 b20 +3c=0

# Plus B2 cond: B2 = lam A2 + d L^2, L=x => L^2=x^2
# B2 = b20 x^2 -2a20 xy -a11/2 y^2
# lam A2 + d x^2 = lam a20 x^2+ lam a11 xy+ lam a02 y^2 + d x^2
# Equate:
# x^2: b20 = lam a20 + d
# xy: -2a20 = lam a11
# y^2: -a11/2 = lam a02

# Collision: A(1,0)= a20 + c = -1, B(1,0)= b20+ cp =0

# So we have unknowns a20,a11,a02,b20,c,cp,lam,d
# Equations:
# (1) -a11^2+4 a02 a20=0
# (2) -2a20 a11 -4 a02 b20=0
# (3) -4a20^2 -2a11 b20 +3c=0
# (4) b20 = lam a20 + d
# (5) -2a20 = lam a11
# (6) -a11/2 = lam a02
# (7) a20 + c = -1
# (8) b20+ cp =0
# That's 8 equations, 8 unknowns (a20,a11,a02,b20,c,cp,lam,d) => potentially solvable!

# Let's solve this system with sympy.

print("\n--- Solving L=x family d=3 ---")
a20_,a11_,a02_,b20_,c_,cp_,lam_,d_ = sp.symbols('a20_ a11_ a02_ b20_ c_ cp_ lam_ d_')
eq1 = -a11_**2 +4*a02_*a20_
eq2 = -2*a20_*a11_ -4*a02_*b20_
eq3 = -4*a20_**2 -2*a11_*b20_ +3*c_
eq4 = b20_ - (lam_*a20_ + d_)
eq5 = -2*a20_ - lam_*a11_
eq6 = -a11_/2 - lam_*a02_
eq7 = a20_ + c_ +1
eq8 = b20_ + cp_

# Solve
sols = sp.solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8], [a20_,a11_,a02_,b20_,c_,cp_,lam_,d_], dict=True)
print(f"Found {len(sols)} solutions for L=x")
for s in sols:
    print(s)

# If no solutions, try solving with Groebner or eliminating
# Let's also attempt general L = x + t y (u=1,v=t) to get more freedom

print("\n--- Attempt general L = x + t y ---")
t = sp.symbols('t')
# Variables: a20,a11,a02,b20,c,cp,lam,d,t
# But t is parameter for L. Let's set u=1, v=t
u_val=1
v_val=t
# Then A2,B2 cond: B2 = lam A2 + d L^2
# L^2 = x^2+2 t xy + t^2 y^2
# Equations:
# b20 = lam a20 + d
# -2a20 = lam a11 + 2 d t
# -a11/2 = lam a02 + d t^2

# J2: (A3)_x+(B3)_y + Jac =0
# A3 = c L^3 = c (x^3+3 t x^2 y+3 t^2 x y^2+ t^3 y^3)
# B3 = cp L^3
# (A3)_x = c *3 L^2 *1 =3c L^2
# (B3)_y = cp *3 L^2 * t =3 cp t L^2
# So (A3)_x+(B3)_y =3(c+ cp t) L^2 =3(c+cp t)(x^2+2 t xy+ t^2 y^2)
# So J2 = 3(c+cp t)L^2 + Jac(A2,B2)=0
# Jac we have expression: (-4a20^2-2 a11 b20) x^2 + (-2a20 a11-4 a02 b20) xy + (-a11^2+4 a02 a20) y^2 = -3(c+cp t)L^2

# So equate coefficients:
# x^2: -4a20^2-2 a11 b20 = -3(c+cp t)
# xy: -2a20 a11-4 a02 b20 = -6 t (c+cp t)
# y^2: -a11^2+4 a02 a20 = -3 t^2 (c+cp t)

# Collision: A(1,0)= a20 + c = -1, B(1,0)= b20+ cp=0

# So we have 3 (B2 cond) +3 (J2) +2 (coll)=8 equations, unknowns a20,a11,a02,b20,c,cp,lam,d,t =9 unknowns => one degree freedom.

# Let's attempt solve with sympy for generic t? Might be heavy. Let's try elimination.

# We'll solve system of 8 equations for 8 unknowns treating t as parameter, then see if solvable.

# Use symbols
a20_,a11_,a02_,b20_,c_,cp_,lam_,d_,t_ = sp.symbols('a20_ a11_ a02_ b20_ c_ cp_ lam_ d_ t_')
eqs = []
eqs.append(b20_ - (lam_*a20_ + d_))  # B2 x2
eqs.append(-2*a20_ - (lam_*a11_ + 2*d_*t_))  # xy
eqs.append(-a11_/2 - (lam_*a02_ + d_*t_**2))  # y2
eqs.append(-4*a20_**2 -2*a11_*b20_ +3*(c_+cp_*t_))  # J2 x2
eqs.append(-2*a20_*a11_ -4*a02_*b20_ +6*t_*(c_+cp_*t_))  # xy
eqs.append(-a11_**2 +4*a02_*a20_ +3*t_**2*(c_+cp_*t_))  # y2
eqs.append(a20_ + c_ +1)
eqs.append(b20_ + cp_)

print("Solving general L family...")
try:
    sols2 = sp.solve(eqs, [a20_,a11_,a02_,b20_,c_,cp_,lam_,d_], dict=True)
    print(f"Found {len(sols2)} parametric solutions in t")
    for s in sols2[:5]:
        print(s)
except Exception as e:
    print("Solve failed:", e)
    # Try Groebner
    try:
        G = sp.groebner(eqs, a20_,a11_,a02_,b20_,c_,cp_,lam_,d_, t_, order='lex')
        print("Groebner:", G)
    except Exception as e2:
        print("Groebner also failed:", e2)

print("\n--- End d=3 analysis ---")
