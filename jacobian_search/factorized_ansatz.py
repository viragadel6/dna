import sympy as sp
x,y=sp.symbols('x y')
v0,w0,s1,s2,t1,t2=sp.symbols('v0 w0 s1 s2 t1 t2')

A = -(1+v0)*x**2 + v0*x**3 + s1*x*y + s2*y**2
B = -w0*x**2 + w0*x**3 + t1*x*y + t2*y**2

# Check A,B in (x,y)^2: yes, no linear, no constant
# Collision: A(1,0)= -(1+v0)+v0 = -1 ok, B(1,0)= -w0 + w0=0 ok

Ax=sp.diff(A,x)
Ay=sp.diff(A,y)
Bx=sp.diff(B,x)
By=sp.diff(B,y)
J=sp.expand(Ax+By+Ax*By-Ay*Bx)
print("J:", J)
poly=sp.Poly(J, x,y)
print(poly.as_dict())
eqs=list(poly.as_dict().values())
print("eqs:", eqs)
# Solve for v0,w0,s1,s2,t1,t2
sols=sp.solve(eqs, [v0,w0,s1,s2,t1,t2], dict=True)
print("Solutions:", sols)

# Try more general: S = s1 x + s2 y + s3 x^2 + s4 xy + s5 y^2 etc for d=4
# Let's try degree 3 for S,T: S = s1 x + s2 y + s3 x^2 + s4 xy + s5 y^2, T = t1 x + t2 y + t3 x^2 + t4 xy + t5 y^2
# Then A = -x^2 + v0 x^2(x-1) + v1 x^2(x-1)*x? Actually v(x) degree1: v(x)=v0+v1 x
# Let's param v(x)=v0+v1 x, w(x)=w0+w1 x

v1,w1,s1,s2,s3,s4,s5,t1,t2,t3,t4,t5=sp.symbols('v1 w1 s1 s2 s3 s4 s5 t1 t2 t3 t4 t5')
v0,w0=sp.symbols('v0 w0')
# A(x,0)= -x^2 + x^2(x-1)(v0+v1 x) = -x^2 + (v0+v1 x)(x^3 - x^2) = -x^2 + v0 x^3 - v0 x^2 + v1 x^4 - v1 x^3 = -(1+v0) x^2 + (v0 - v1) x^3 + v1 x^4
# B(x,0)= x^2(x-1)(w0+w1 x)= w0 x^3 - w0 x^2 + w1 x^4 - w1 x^3 = -w0 x^2 + (w0 - w1) x^3 + w1 x^4

# S = s1 x + s2 y + s3 x^2 + s4 x y + s5 y^2
# T = t1 x + t2 y + t3 x^2 + t4 x y + t5 y^2

S = s1*x + s2*y + s3*x**2 + s4*x*y + s5*y**2
T = t1*x + t2*y + t3*x**2 + t4*x*y + t5*y**2

A2 = -(1+v0)*x**2 + (v0 - v1)*x**3 + v1*x**4 + y*S
B2 = -w0*x**2 + (w0 - w1)*x**3 + w1*x**4 + y*T

# Expand
A = sp.expand(A2)
B = sp.expand(B2)
print("\n\nDegree 4 factorized ansatz")
print("A:", A)
print("B:", B)
Ax=sp.diff(A,x)
Ay=sp.diff(A,y)
Bx=sp.diff(B,x)
By=sp.diff(B,y)
J=sp.expand(Ax+By+Ax*By-Ay*Bx)
poly=sp.Poly(J, x,y)
print(f"J degree: {poly.total_degree()}, terms: {len(poly.as_dict())}")
eqs=list(poly.as_dict().values())
print(f"Number eqs: {len(eqs)}")
# Variables: v0,v1,w0,w1,s1..s5,t1..t5 = 2+2+5+5=14 vars, eqs maybe ~? Let's count monomials degree up to 6? Since A,B degree4, J degree up to 6, number monomials degree1..6 =27, so 27 eqs, 14 vars overdetermined
# Try to solve with Groebner? Might be heavy
# Try to find solution with some vars zero to reduce

# Let's try to set many s,t zero: s1=s2=t1=t2=0? Actually S must have no constant, but we can set linear part zero to keep A,B in (x,y)^2? Wait y*S with S linear gives quadratic terms, okay. If we set s1=s2=0, then y*S = y*(quadratic) = degree3, so A would have no quadratic terms in y? But A must have degree at least2, but -x^2 already provides, okay.
# Let's try sparse: S = s3 x^2, T = t3 x^2

s1_val=0
s2_val=0
s4_val=0
s5_val=0
t1_val=0
t2_val=0
t4_val=0
t5_val=0

S_sparse = s3*x**2
T_sparse = t3*x**2

A_sparse = -(1+v0)*x**2 + (v0 - v1)*x**3 + v1*x**4 + y*S_sparse
B_sparse = -w0*x**2 + (w0 - w1)*x**3 + w1*x**4 + y*T_sparse

Ax=sp.diff(A_sparse,x)
Ay=sp.diff(A_sparse,y)
Bx=sp.diff(B_sparse,x)
By=sp.diff(B_sparse,y)
J_sparse=sp.expand(Ax+By+Ax*By-Ay*Bx)
poly_sparse=sp.Poly(J_sparse, x,y)
print("\nSparse J:", poly_sparse.as_dict())
eqs_sparse=list(poly_sparse.as_dict().values())
print("eqs_sparse:", eqs_sparse)
# Solve for v0,v1,w0,w1,s3,t3
sols_sparse=sp.solve(eqs_sparse, [v0,v1,w0,w1,s3,t3], dict=True)
print("Solutions sparse:", sols_sparse)

# Try even sparser: A = -x^2 + v0 x^3 + y * s1 x, B = w0 x^3 + y * t1 x? Let's brute force search for low-degree solutions with random support

print("\nDone")
