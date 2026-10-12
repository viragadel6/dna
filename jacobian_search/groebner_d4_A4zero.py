import sympy as sp

x,y=sp.symbols('x y')

# Variables for A4=0 case as before
a20,a11,a02,b20, a30,a21,a12,a03,b30, b40 = sp.symbols('a20 a11 a02 b20 a30 a21 a12 a03 b30 b40')

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

B4 = b40*x**4 + b31_expr*x**3*y + b22_expr*x**2*y**2 + b13_expr*x*y**3 + b04_expr*y**4

JacA2B4 = sp.diff(A2,x)*sp.diff(B4,y) - sp.diff(A2,y)*sp.diff(B4,x)
JacA3B3 = sp.diff(A3,x)*sp.diff(B3,y) - sp.diff(A3,y)*sp.diff(B3,x)
J4 = sp.expand(JacA2B4 + JacA3B3)
polyJ4 = sp.Poly(J4, x,y)
eqsJ4 = list(polyJ4.as_dict().values())

JacA3B4 = sp.diff(A3,x)*sp.diff(B4,y) - sp.diff(A3,y)*sp.diff(B4,x)
J5 = sp.expand(JacA3B4)
polyJ5 = sp.Poly(J5, x,y)
eqsJ5 = list(polyJ5.as_dict().values())

A_coll = (A2 + A3).subs({x:1,y:0})
B_coll = (B2 + B3 + B4).subs({x:1,y:0})
eq_coll1 = A_coll + 1
eq_coll2 = B_coll

all_eqs = eqsJ4 + eqsJ5 + [eq_coll1, eq_coll2]
print(f"Total eqs: {len(all_eqs)}")
print("Vars: a20,a11,a02,b20,a30,a21,a12,a03,b30,b40")

# Try Groebner with lex order, but may be heavy
# Let's try with groebner and see if 1 is in basis (inconsistent)
try:
    G = sp.groebner(all_eqs, a20,a11,a02,b20,a30,a21,a12,a03,b30,b40, order='lex')
    print("Groebner basis computed")
    print(G)
    # Check if 1 in basis
    if sp.Poly(1, a20,a11,a02,b20,a30,a21,a12,a03,b30,b40) in G.polys:
        print("System inconsistent (contains 1)")
    else:
        # Check if basis contains constant
        for poly in G.polys:
            if poly.is_number and poly!=0:
                print("Constant non-zero in basis => inconsistent")
                break
        print("Groebner does not contain 1, may have solutions")
except Exception as e:
    print(f"Groebner failed: {e}")

# Try with grlex order (maybe faster)
try:
    G2 = sp.groebner(all_eqs, a20,a11,a02,b20,a30,a21,a12,a03,b30,b40, order='grlex')
    print("\nGroebner grlex computed")
    print(G2)
except Exception as e:
    print(f"Groebner grlex failed: {e}")

print("\nDone")
