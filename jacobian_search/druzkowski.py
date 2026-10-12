import sympy as sp
a,b,c,d=sp.symbols('a b c d')
eq1 = a**3 + d*c**2
eq2 = a**2*b + c*d**2
eq3 = a*b**2 + d**3
eq4 = a*d - b*c
sols=sp.solve([eq1,eq2,eq3,eq4],[a,b,c,d],dict=True)
print(sols)
# Solve more generally allowing parameters
# Use groebner
G=sp.groebner([eq1,eq2,eq3,eq4],a,b,c,d,order='lex')
print(G)
# Let's find all solutions with some free params
# From eq4, if a!=0, d= b c /a
# Plug into eq1: a**3 + (b c /a) c**2 = a**3 + b c**3 /a =0 => a**4 + b c**3=0
# eq2: a**2 b + c (b c /a)**2 = a**2 b + c b**2 c**2 /a**2 = a**2 b + b**2 c**3 /a**2=0 => multiply a**2: a**4 b + b**2 c**3= b(a**4+ b c**3)=0
# eq3: a b**2 + (b c /a)**3 = a b**2 + b**3 c**3 /a**3 =0 => multiply a**3: a**4 b**2 + b**3 c**3 = b**2(a**4+ b c**3)=0
# So all reduce to a**4 + b c**3=0
# So solutions: choose b,c arbitrary, set a = (-b c**3)**(1/4), d= b c /a
# Then we get non-zero a,d
# Let's pick b=1, c=1 => a^4 = -1 => a = exp(i pi/4) etc complex allowed! Since we are over C, we can have complex solutions!
# So over C, there are non-trivial Druzkowski maps with both components non-zero!
# For example, b=1,c=1 => a^4 = -1 => a = (1+i)/sqrt(2) ??? Actually (-1)^(1/4) = (1+i)/sqrt(2) *? Let's compute: (-1)^(1/4) = e^{i pi/4}= (1+i)/sqrt2, e^{i 3pi/4}=(-1+i)/sqrt2, etc.
# Then d = b c /a =1/a
# So we have a map: F1 = x + (a x + y)^3, F2 = y + (x + d y)^3 with a^4=-1, d=1/a
# Does this satisfy trace zero? Let's test numeric

import cmath
b_val=1
c_val=1
a_val = complex(0,1)**0.5 # sqrt(i) = e^{i pi/4}
# Actually (-1)**0.25 = e^{i pi/4}
a_val = cmath.exp(1j*cmath.pi/4)
d_val = b_val*c_val/a_val
print(f"a={a_val}, d={d_val}, b={b_val}, c={c_val}")
# Check eqs
print(a_val**3 + d_val*c_val**2)
print(a_val**2*b_val + c_val*d_val**2)
print(a_val*b_val**2 + d_val**3)
print(a_val*d_val - b_val*c_val)

# So over C, non-trivial Druzkowski maps exist! Over R they don't (since a^4 + b c^3 =0 with real b,c => a real => a^4 >=0, b c^3 could be negative, possible? If b c^3 negative, a^4 = -b c^3 positive, so real solutions exist too! For b=1,c=-1 => b c^3 = -1 => a^4=1 => a=1 => d= b c /a = -1 => check eq1: 1 + (-1)*1=0 ok, eq2:1 + (-1)*1=0, eq3:1 + (-1)=0. So a=1,b=1,c=-1,d=-1 gives real solution! Let's test
a_val=1
b_val=1
c_val=-1
d_val=-1
print("\nReal test a=1,b=1,c=-1,d=-1")
print(a_val**3 + d_val*c_val**2)
print(a_val**2*b_val + c_val*d_val**2)
print(a_val*b_val**2 + d_val**3)
print(a_val*d_val - b_val*c_val)

# So F1 = x + (x + y)^3, F2 = y + (-x - y)^3 = y - (x+y)^3
# Then F = (x + (x+y)^3, y - (x+y)^3)
# Jacobian: compute
x,y=sp.symbols('x y')
F1 = x + (x+y)**3
F2 = y - (x+y)**3
J = sp.Matrix([[sp.diff(F1,x), sp.diff(F1,y)],[sp.diff(F2,x), sp.diff(F2,y)]]).det()
print("Det for real example:", sp.expand(J))
# This is Keller! Does it have collision? Need F(p)=F(q)
# F1 = x + s^3, F2 = y - s^3 where s=x+y
# So F1+F2 = x+y = s
# So s = F1+F2
# Then F1 = x + s^3 => x = F1 - s^3 = F1 - (F1+F2)^3
# Similarly y = F2 + (F1+F2)^3
# So inverse exists polynomial! So it's automorphism, injective.

# In general, Druzkowski maps with JH nilpotent are known to be invertible? The conjecture is that all such are invertible, but not proven in general.

# So this family won't give counterexample.

# Let's search for more general maps where H homogeneous degree 3 but not necessarily of form (linear)^3, but general homogeneous degree 3 with nilpotent Jacobian.

# General homogeneous degree 3 H = (H1,H2) with H1,H2 homogeneous degree3, JH nilpotent: trace=0 and det=0.

# Trace = (H1)_x + (H2)_y =0, det = (H1)_x (H2)_y - (H1)_y (H2)_x =0

# We can param general H1 = sum a_i x^{3-i} y^i, H2 = sum b_i x^{3-i} y^i
# Then (H1)_x = 3a0 x^2+2a1 xy+ a2 y^2, (H2)_y = b1 x^2+2b2 xy+3b3 y^2
# Trace zero: 3a0 + b1 =0, 2a1+2b2=0 => a1+ b2=0, a2+3b3=0
# So b1=-3a0, b2=-a1, b3= -a2/3

# So H2 = b0 x^3 -3a0 x^2 y -a1 x y^2 -a2/3 y^3

# Then det condition: (H1)_x (H2)_y - (H1)_y (H2)_x =0

# Compute with sympy general

a0,a1,a2,a3,b0 = sp.symbols('a0 a1 a2 a3 b0')
# H1 = a0 x^3 + a1 x^2 y + a2 x y^2 + a3 y^3
# H2 = b0 x^3 + b1 x^2 y + b2 x y^2 + b3 y^3 with b1,b2,b3 as above
b1 = -3*a0
b2 = -a1
b3 = -a2/3
H1 = a0*x**3 + a1*x**2*y + a2*x*y**2 + a3*y**3
H2 = b0*x**3 + b1*x**2*y + b2*x*y**2 + b3*y**3
H1x = sp.diff(H1,x)
H1y = sp.diff(H1,y)
H2x = sp.diff(H2,x)
H2y = sp.diff(H2,y)
det = sp.expand(H1x*H2y - H1y*H2x)
print("\ndet general homogeneous degree3 nilpotent:")
print(det)
poly_det = sp.Poly(det, x,y)
print(poly_det.as_dict())
# Solve det=0 identically => all coefficients zero
eqs_det = list(poly_det.as_dict().values())
print("eqs_det:", eqs_det)
sols_det = sp.solve(eqs_det, [a0,a1,a2,a3,b0], dict=True)
print("Solutions for homogeneous degree3 Keller:", sols_det)
# This is system of 5 equations (degree4 homogeneous has 5 coeffs) in 5 unknowns
# Solve

# Let's try groebner
try:
    G2 = sp.groebner(eqs_det, a0,a1,a2,a3,b0, order='lex')
    print("Groebner:", G2)
except Exception as e:
    print("Groebner failed:", e)

# Let's attempt to solve with some free param: maybe set a3=0?
# Try numeric solving: use elimination

# Let's attempt to solve generally: from det expression, we can try to factor
# det = ?

# Let's attempt to find all solutions where H1 = c L^3, H2 = d L^3 as earlier, does that satisfy trace zero? Let's see

# Suppose L = u x + v y, H1 = c L^3, H2 = d L^3
# Then H1x =3c u L^2, H1y=3c v L^2, H2x=3d u L^2, H2y=3d v L^2
# Trace =3c u L^2 +3d v L^2 =3L^2(c u + d v)=0 => c u + d v=0
# Det = (3c u L^2)(3d v L^2)-(3c v L^2)(3d u L^2)=0 automatically
# So condition is c u + d v=0, which matches earlier linear form analysis (p = c L^3, q=d L^3, u p + v q? Actually p' =3c L^2, q'=3d L^2, u p' + v q' =3L^2(c u + d v)=0)
# So this family is included

# Are there other solutions where H1,H2 not proportional L^3? Let's try to solve det equations

# Use sympy solve with some assumptions: maybe set b0=0 to simplify
eqs_det_b0_0 = [e.subs({b0:0}) for e in eqs_det]
print("\nDet eqs with b0=0:", eqs_det_b0_0)
sols_b0_0 = sp.solve(eqs_det_b0_0, [a0,a1,a2,a3], dict=True)
print("Sols b0=0:", sols_b0_0)

# Let's brute force search for solutions with random integer values for a0,a1,a2,a3,b0 small and check det zero?
# Could attempt to find parametrization: det=0 means H1,H2 algebraically dependent, so they share common factor L? Since degree3, if they are dependent, they must be of form c L^3 and d L^3? Actually earlier we argued for same degree homogeneous, Jac=0 => proportional. So for homogeneous degree3, Jac=0 => H1 and H2 proportional! Wait we earlier proved for binary homogeneous same degree, Jac=0 => proportional. So det=0 (which is Jac(H1,H2)=0) implies H1 and H2 proportional! Therefore H2 = k H1 for some k.

# Let's check: If H2 = k H1, then trace = (H1)_x + k (H1)_y =0? Actually (H2)_y = k (H1)_y, so trace = (H1)_x + k (H1)_y =0. This is a PDE for H1.

# So we can set H2 = k H1, then condition is (H1)_x + k (H1)_y =0 => H1 is function of (x - k y)? Let's solve: PDE f_x + k f_y =0 => f = g(y - k x?) Actually characteristic: dx/1 = dy/k => dy = k dx => y - k x = const, so f = h(y - k x) ??? Let's check: f(x,y)= h(y - k x) => f_x = -k h', f_y = h' => f_x + k f_y = -k h' + k h' =0. Yes.

# So H1 = h(y - k x) homogeneous degree3 => h(t) = c t^3, and argument y - k x is linear, so H1 = c (y - k x)^3, H2 = k c (y - k x)^3

# Then F1 = x + c (y - k x)^3, F2 = y + k c (y - k x)^3
# This is of the form earlier with L = y - k x, and c u + d v =0? Let's see: L = -k x + y, u=-k, v=1, c1=c, d1=k c => c1 u + d1 v = c(-k)+ k c*1=0 yes.

# So the only homogeneous degree3 Keller maps are of this L^3 proportional form! Therefore they are automorphisms (since triangular after linear change).

# Thus degree 3 Keller maps are all automorphisms, no counterexample.

# This matches known result: homogeneous degree 3 Keller in 2D are tame.

# So we need degree >=4.

print("\nConclusion: degree 3 homogeneous Keller maps are exactly L^3 family, no counterexample")
