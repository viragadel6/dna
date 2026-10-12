import sympy as sp

a20,a11,a02,c,k,d,v = sp.symbols('a20 a11 a02 c k d v')
# u=1
# Equations E1..E6 as derived
E1 = 2*a20 + k*a11 + 2*d*v
E2 = a11 + 2*k*a02 + 2*d*v**2
E3 = 3*c*(1+ k*v) + 2*d*(2*a20*v - a11)
E4 = 3*c*(1+ k*v)*v + 2*d*(a11*v - 2*a02)
E5 = a20 + c + 1
E6 = k*a20 + d + k*c

print("Solving case3 with u=1")
# From E5, E6 => d=k, c = -1-a20
# Let's solve generally using sympy solve
sols = sp.solve([E1,E2,E3,E4,E5,E6], [a20,a11,a02,c,d], dict=True)
print(f"Solutions parametric in k,v: {len(sols)}")
for s in sols:
    print(s)

# The above solve eliminated k,v? Let's solve for a20,a11,a02,c,d in terms of k,v
# Use solve with k,v as parameters: solve linear system E1,E2,E5,E6 for a20,a11,a02,c,d? E1,E2 linear in a11,a02, E5,E6 linear in c,d.
# Let's do manual elimination as derived: d=k, c=-1-a20
# Then E1 => a11 = (-2a20 -2k v)/k
# E2 => a02 = -(a11+2k v^2)/(2k)
# Plug into E3,E4 to get equations for a20,k,v

a20_sym,k_sym,v_sym = sp.symbols('a20_sym k_sym v_sym')
c_expr = -1 - a20_sym
d_expr = k_sym
a11_expr = (-2*a20_sym -2*k_sym*v_sym)/k_sym
a02_expr = -(a11_expr + 2*k_sym*v_sym**2)/(2*k_sym)

E3_sub = 3*c_expr*(1+ k_sym*v_sym) + 2*d_expr*(2*a20_sym*v_sym - a11_expr)
E4_sub = 3*c_expr*(1+ k_sym*v_sym)*v_sym + 2*d_expr*(a11_expr*v_sym - 2*a02_expr)

print("\nE3_sub:", sp.expand(E3_sub))
print("E4_sub:", sp.expand(E4_sub))

# Simplify E3_sub, E4_sub
E3_s = sp.simplify(E3_sub)
E4_s = sp.simplify(E4_sub)
print("E3_s simplified:", E3_s)
print("E4_s simplified:", E4_s)

# Solve E3_s=0, E4_s=0 for a20 in terms of k,v
sol_a20_E3 = sp.solve(E3_s, a20_sym)
print("\nSolve E3 for a20:", sol_a20_E3)
for sol in sol_a20_E3:
    E4_check = E4_s.subs({a20_sym:sol})
    print(f"a20={sol}, E4 becomes {sp.simplify(E4_check)}")
    # Solve E4_check=0 for k,v?
    # E4_check is expression in k,v
    # Solve for v?
    # Let's factor
    print("Factored E4_check:", sp.factor(E4_check))

# Let's solve system E3_s=0, E4_s=0 simultaneously for a20 and maybe k or v
sols2 = sp.solve([E3_s, E4_s], [a20_sym], dict=True)
print("\nSolutions for a20 from both eqs:", sols2)

# Solve for a20,k maybe
# Use groebner for E3_s, E4_s in a20,k with v param? Let's treat v as param, solve for a20,k
# From E3_s: expression linear in a20? Let's check degree
print("\nDegree of E3_s in a20:", sp.Poly(E3_s, a20_sym).degree())
print("Degree of E4_s in a20:", sp.Poly(E4_s, a20_sym).degree())

# They are linear in a20, so we can solve both and equate
# Solve E3_s=0 => a20 = ...
a20_from_E3 = sp.solve(E3_s, a20_sym)[0]
print("a20 from E3:", a20_from_E3)
a20_from_E4 = sp.solve(E4_s, a20_sym)[0]
print("a20 from E4:", a20_from_E4)
# Set equal
eq_equal = sp.simplify(a20_from_E3 - a20_from_E4)
print("Equal condition (a20_E3 - a20_E4)=0:", sp.simplify(eq_equal))
print("Factored:", sp.factor(eq_equal))
# Solve eq_equal=0 for k,v
# Multiply by denominators
num, den = sp.together(eq_equal).as_numer_denom()
print("Numerator:", sp.factor(num))
# Solve numerator=0 for k
sol_k = sp.solve(num, k_sym)
print("Solve for k:", sol_k)
# For each k, find a20
for ks in sol_k:
    a20_val = a20_from_E3.subs({k_sym:ks})
    print(f"k={ks}, a20={sp.simplify(a20_val)}")

# Let's also solve numerator=0 for v
# The numerator factor might give relation between k and v
# Let's factor numerator
# Try to solve for v in terms of k
# Use sympy solve for v
try:
    sol_v = sp.solve(num, v_sym)
    print("Solve for v:", sol_v)
except Exception as e:
    print("Solve for v failed:", e)

# Let's attempt to find any solution with k!=0
# Brute force: pick k=1, solve for v and a20
k_val=1
E3_k1 = E3_s.subs({k_sym:k_val})
E4_k1 = E4_s.subs({k_sym:k_val})
print(f"\nFor k=1: E3={E3_k1}, E4={E4_k1}")
sol_a20_k1_E3 = sp.solve(E3_k1, a20_sym)
print("a20 from E3 k=1:", sol_a20_k1_E3)
for a20v in sol_a20_k1_E3:
    print("Check E4:", E4_k1.subs({a20_sym:a20v}), "factor:", sp.factor(E4_k1.subs({a20_sym:a20v})))
    # Solve for v
    eq_v = E4_k1.subs({a20_sym:a20v})
    sol_v_k1 = sp.solve(eq_v, v_sym)
    print("sol v:", sol_v_k1)

# Try k=0 case separately
print("\n--- k=0 case ---")
# k=0 => B3 =0, case2
# Then our earlier derivation with d=k=0 => d=0
# Equations:
# E1: 2a20=0 => a20=0
# E2: a11=0
# E3: 3c*1 +0=0 => c=0
# But c = -1 -a20 = -1 => contradiction, so k=0 no solution in this branch (case3). So case3 requires k!=0.
# Let's analyze case2: B3=0, A3=c L^3, B2 = d L^2 (since k=0, B2 = d L^2)
# J1: 2a20 +2d u v=0, a11+2d v^2=0
# J2: 3c u L? Actually (A3)_x =3c u L^2, Jac(A2,B2)=d Jac(A2,L^2)=2d L Jac(A2,L)
# So J2: 3c u L^2 +2d L Jac(A2,L)=0 => L[3c u L +2d Jac(A2,L)]=0 => 3c u L +2d Jac=0
# Similar as before with k=0 but with u factor

# Let's do case2 separately

print("\n--- Case2: B3=0, A3=c L^3, B2=d L^2 ---")
# Variables a20,a11,a02,c,d,v with u=1
a20_,a11_,a02_,c_,d_,v_ = sp.symbols('a20_ a11_ a02_ c_ d_ v_')
# J1: B2 = d L^2 = d(x^2+2v xy+v^2 y^2) => b20=d, b11=2d v, b02=d v^2
# J1 conditions: 2a20 + b11=0 => 2a20+2d v=0
# a11+2b02=0 => a11+2d v^2=0
E1_c2 = 2*a20_ + 2*d_*v_
E2_c2 = a11_ + 2*d_*v_**2
# J2 bracket: 3c u L +2d Jac(A2,L)=0, u=1, L=x+v y, Jac(A2,L)= x(2a20 v -a11)+ y(a11 v -2a02)
# So: 3c (x+v y) +2d[ x(2a20 v -a11)+ y(a11 v -2a02)]=0
# Coeff x: 3c +2d(2a20 v -a11)=0
# Coeff y: 3c v +2d(a11 v -2a02)=0
E3_c2 = 3*c_ + 2*d_*(2*a20_*v_ - a11_)
E4_c2 = 3*c_*v_ + 2*d_*(a11_*v_ -2*a02_)
# Collision: A(1,0)= a20 +c = -1, B(1,0)= b20 = d =0? Wait B = B2 = d L^2, so B(1,0)=d*1= d must be 0 for collision B(1,0)=0
E5_c2 = a20_ + c_ +1
E6_c2 = d_  # =0

print("Case2 eqs:", E1_c2, E2_c2, E3_c2, E4_c2, E5_c2, E6_c2)
sols_c2 = sp.solve([E1_c2,E2_c2,E3_c2,E4_c2,E5_c2,E6_c2], [a20_,a11_,a02_,c_,d_,v_], dict=True)
print("Solutions case2:", sols_c2)
# From E6, d=0 => then E1=>a20=0, E2=>a11=0, E3=>3c=0=>c=0, E5=>0+0+1=0 contradiction => no solution

print("\n--- Case1: A3=0, B3=cp L^3, A2 =? Actually A3=0, B3=cp L^3, Jac(A2,B3)=0 => A2 = a L^2 ---")
# A3=0, B3=cp L^3, A2 = a L^2? Let's derive
# A3=0, B3=cp L^3, condition Jac(A2,B3)=0 => Jac(A2,L)=0 => A2 = a L^2
# So A2 = a L^2, B2 =? B2 general but J1 and J2 conditions
# J1: A2 = a L^2 = a(x^2+2v xy+v^2 y^2) => a20=a, a11=2a v, a02=a v^2
# B2 = b20 x^2 + b11 xy + b02 y^2 with J1: 2a20 + b11=0 => 2a + b11=0 => b11=-2a
# a11+2b02=0 =>2a v +2b02=0 => b02= -a v
# So B2 = b20 x^2 -2a xy -a v y^2
# J2: (B3)_y + Jac(A2,B2)=0
# (B3)_y = 3cp v L^2
# Jac(A2,B2)=? Since A2 = a L^2, Jac(A2,B2)= a Jac(L^2,B2)=2a L Jac(L,B2)
# Compute Jac(L,B2): L=x+v y, B2_x=2b20 x -2a y, B2_y= -2a x -2a v y? Wait B2_y = b11 x +2b02 y = -2a x -2a v y
# So Jac(L,B2)= L_x B2_y - L_y B2_x =1*B2_y - v*B2_x = (-2a x -2a v y) - v(2b20 x -2a y)= -2a x -2a v y -2b20 v x +2a v y = -2a x -2b20 v x = -2x(a+ b20 v)
# So Jac(A2,B2)=2a L * (-2x(a+ b20 v)) = -4a x L (a+ b20 v)
# This is degree3? Wait L is degree1, x is degree1 => degree2? Actually L * x = degree2, so Jac is degree2? But earlier we expected degree2, okay.
# J2 = 3cp v L^2 + Jac(A2,B2)=0 => 3cp v L^2 -4a x L (a+ b20 v)=0 => L[3cp v L -4a x(a+ b20 v)]=0 => 3cp v L -4a x(a+ b20 v)=0
# So: 3cp v (x+v y) -4a x(a+ b20 v)=0 => x[3cp v -4a(a+ b20 v)] + y[3cp v^2]=0
# So coefficients:
# y: 3cp v^2=0 => either cp=0 or v=0
# x: 3cp v -4a(a+ b20 v)=0
# Collision: A(1,0)= A2(1,0)=a = -1, B(1,0)= B2(1,0)+B3(1,0)= b20 + cp =0
# So a=-1
# Then conditions:
# If v=0: then y coeff automatically 0, x coeff: 0 -4(-1)(-1+ b20*0)= -4*1= -4? Wait compute: a=-1, so -4a(a+ b20 v)= -4(-1)(-1+0)= -4*1= -4 !=0 => no solution
# If cp=0: then B3=0, both A3,B3 zero => degree2 case, which we know no collision? Let's see: cp=0 => B3=0, then J2 becomes Jac(A2,B2)=0, but we already have A2=a L^2, B2 as above, Jac(A2,B2)= -4a x L(a+ b20 v)=0 => requires a+ b20 v=0 => -1+ b20 v=0 => b20=1/v if v!=0
# Collision: b20=0 (since cp=0, B(1,0)=b20=0) => contradiction with b20=1/v
# So case1 no solution

# Let's verify with sympy for case1
a_,b20_,cp_,v_ = sp.symbols('a_ b20_ cp_ v_')
# a = A2 coeff, a20=a, a11=2a v, a02=a v^2
# B2: b20, b11=-2a, b02=-a v
# J2 conditions derived above but let's compute directly with sympy for general case1
x,y=sp.symbols('x y')
a_sym,b20_sym,cp_sym,v_sym=sp.symbols('a_sym b20_sym cp_sym v_sym')
L_sym = x + v_sym*y
A2_sym = a_sym*L_sym**2
B2_sym = b20_sym*x**2 -2*a_sym*x*y - a_sym*v_sym*y**2
B3_sym = cp_sym*L_sym**3
# J1 already satisfied by construction? Check J1: (A2)_x+(B2)_y
J1_sym = sp.diff(A2_sym,x)+sp.diff(B2_sym,y)
print("\nCase1 J1:", sp.expand(J1_sym))  # Should be zero? Let's see
# Compute J2: (B3)_y + Jac(A2,B2)
Jac_A2_B2_sym = sp.diff(A2_sym,x)*sp.diff(B2_sym,y)-sp.diff(A2_sym,y)*sp.diff(B2_sym,x)
J2_sym = sp.expand(sp.diff(B3_sym,y)+Jac_A2_B2_sym)
print("Case1 J2:", J2_sym)
polyJ2 = sp.Poly(J2_sym, x,y)
print("J2 coeffs:", polyJ2.as_dict())
# Collision: A(1,0)=a=-1, B(1,0)=b20+cp=0
# So set a=-1
J2_a_1 = J2_sym.subs({a_sym:-1})
print("J2 with a=-1:", sp.expand(J2_a_1))
polyJ2_a1 = sp.Poly(J2_a_1, x,y)
print(polyJ2_a1.as_dict())
# Solve J2=0 for b20,cp,v
sols_J2 = sp.solve(list(polyJ2_a1.as_dict().values()), [b20_sym,cp_sym], dict=True)
print("Solve J2 for b20,cp with a=-1:", sols_J2)
# For each solution, check collision b20+cp=0?
# Let's solve including collision
eqs_case1 = list(polyJ2_a1.as_dict().values()) + [b20_sym+cp_sym]
sols_case1_full = sp.solve(eqs_case1, [b20_sym,cp_sym,v_sym], dict=True)
print("Full case1 with collision:", sols_case1_full)

print("\nDone case analysis")
