import sympy as sp, random

x,y=sp.symbols('x y')

# Degree 5: we will generate general A,B up to degree5 with 18 monomials each, but we will eliminate linear J1,J2,J3,J4 as before to get 22 free vars
# Let's do step by step

# Variables for free after elimination:
# A2,B2: a20,a11,a02,b20 (4)
# A3,B3: a30,a21,a12,a03,b30 (5) total 9
# A4,B4: a40,a31,a22,a13,a04,b40 (6) total 15
# A5,B5: a50,a41,a32,a23,a14,a05,b50 (7) total 22

a20,a11,a02,b20, a30,a21,a12,a03,b30, a40,a31,a22,a13,a04,b40, a50,a41,a32,a23,a14,a05,b50 = sp.symbols('a20 a11 a02 b20 a30 a21 a12 a03 b30 a40 a31 a22 a13 a04 b40 a50 a41 a32 a23 a14 a05 b50')

# A2,B2
b11 = -2*a20
b02 = -a11/2
A2 = a20*x**2 + a11*x*y + a02*y**2
B2 = b20*x**2 + b11*x*y + b02*y**2

# A3
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

# A4
A4 = a40*x**4 + a31*x**3*y + a22*x**2*y**2 + a13*x*y**3 + a04*y**4
JacA2B3 = sp.diff(A2,x)*sp.diff(B3,y) - sp.diff(A2,y)*sp.diff(B3,x)
JacA3B2 = sp.diff(A3,x)*sp.diff(B2,y) - sp.diff(A3,y)*sp.diff(B2,x)
J3_rhs = sp.expand(JacA2B3 + JacA3B2)
polyJ3 = sp.Poly(J3_rhs, x,y)
c_x3 = polyJ3.coeff_monomial(x**3)
c_x2y = polyJ3.coeff_monomial(x**2*y)
c_xy2 = polyJ3.coeff_monomial(x*y**2)
c_y3 = polyJ3.coeff_monomial(y**3)
b31_expr = -4*a40 - c_x3
b22_expr = (-3*a31 - c_x2y)/2
b13_expr = (-2*a22 - c_xy2)/3
b04_expr = (-a13 - c_y3)/4
B4 = b40*x**4 + b31_expr*x**3*y + b22_expr*x**2*y**2 + b13_expr*x*y**3 + b04_expr*y**4

# A5
A5 = a50*x**5 + a41*x**4*y + a32*x**3*y**2 + a23*x**2*y**3 + a14*x*y**4 + a05*y**5
JacA2B4 = sp.diff(A2,x)*sp.diff(B4,y) - sp.diff(A2,y)*sp.diff(B4,x)
JacA3B3 = sp.diff(A3,x)*sp.diff(B3,y) - sp.diff(A3,y)*sp.diff(B3,x)
JacA4B2 = sp.diff(A4,x)*sp.diff(B2,y) - sp.diff(A4,y)*sp.diff(B2,x)
J4_rhs = sp.expand(JacA2B4 + JacA3B3 + JacA4B2)
polyJ4 = sp.Poly(J4_rhs, x,y)
# J4: (A5)_x + (B5)_y = -J4_rhs
# A5x = 5a50 x^4 +4a41 x^3 y +3a32 x^2 y^2 +2a23 x y^3 + a14 y^4
# B5y = b41 x^4 +2b32 x^3 y +3b23 x^2 y^2 +4b14 x y^3 +5b05 y^4
# So we can solve for B5 coeffs b41,b32,b23,b14,b05 in terms of A5 and J4_rhs
# Let's get coeffs of J4_rhs degree4
c_x4 = polyJ4.coeff_monomial(x**4)
c_x3y = polyJ4.coeff_monomial(x**3*y)
c_x2y2 = polyJ4.coeff_monomial(x**2*y**2)
c_xy3 = polyJ4.coeff_monomial(x*y**3)
c_y4 = polyJ4.coeff_monomial(y**4)

# Then:
# 5a50 + b41 + c_x4 =0 => b41 = -5a50 - c_x4
# 4a41 +2b32 + c_x3y=0 => b32 = (-4a41 - c_x3y)/2
# 3a32 +3b23 + c_x2y2=0 => b23 = -a32 - c_x2y2/3
# 2a23 +4b14 + c_xy3=0 => b14 = (-2a23 - c_xy3)/4
# a14 +5b05 + c_y4=0 => b05 = (-a14 - c_y4)/5

b41_expr = -5*a50 - c_x4
b32_expr = (-4*a41 - c_x3y)/2
b23_expr = -a32 - c_x2y2/3
b14_expr = (-2*a23 - c_xy3)/4
b05_expr = (-a14 - c_y4)/5

B5 = b50*x**5 + b41_expr*x**4*y + b32_expr*x**3*y**2 + b23_expr*x**2*y**3 + b14_expr*x*y**4 + b05_expr*y**5

# Now remaining equations: J5,J6,J7,J8
JacA2B5 = sp.diff(A2,x)*sp.diff(B5,y) - sp.diff(A2,y)*sp.diff(B5,x)
JacA3B4 = sp.diff(A3,x)*sp.diff(B4,y) - sp.diff(A3,y)*sp.diff(B4,x)
JacA4B3 = sp.diff(A4,x)*sp.diff(B3,y) - sp.diff(A4,y)*sp.diff(B3,x)
JacA5B2 = sp.diff(A5,x)*sp.diff(B2,y) - sp.diff(A5,y)*sp.diff(B2,x)
J5 = sp.expand(JacA2B5 + JacA3B4 + JacA4B3 + JacA5B2)
polyJ5 = sp.Poly(J5, x,y)
eqsJ5 = list(polyJ5.as_dict().values())

JacA3B5 = sp.diff(A3,x)*sp.diff(B5,y) - sp.diff(A3,y)*sp.diff(B5,x)
JacA4B4 = sp.diff(A4,x)*sp.diff(B4,y) - sp.diff(A4,y)*sp.diff(B4,x)
JacA5B3 = sp.diff(A5,x)*sp.diff(B3,y) - sp.diff(A5,y)*sp.diff(B3,x)
J6 = sp.expand(JacA3B5 + JacA4B4 + JacA5B3)
polyJ6 = sp.Poly(J6, x,y)
eqsJ6 = list(polyJ6.as_dict().values())

JacA4B5 = sp.diff(A4,x)*sp.diff(B5,y) - sp.diff(A4,y)*sp.diff(B5,x)
JacA5B4 = sp.diff(A5,x)*sp.diff(B4,y) - sp.diff(A5,y)*sp.diff(B4,x)
J7 = sp.expand(JacA4B5 + JacA5B4)
polyJ7 = sp.Poly(J7, x,y)
eqsJ7 = list(polyJ7.as_dict().values())

JacA5B5 = sp.diff(A5,x)*sp.diff(B5,y) - sp.diff(A5,y)*sp.diff(B5,x)
J8 = sp.expand(JacA5B5)
polyJ8 = sp.Poly(J8, x,y)
eqsJ8 = list(polyJ8.as_dict().values())

A_coll = (A2 + A3 + A4 + A5).subs({x:1,y:0})
B_coll = (B2 + B3 + B4 + B5).subs({x:1,y:0})
eq_coll1 = A_coll + 1
eq_coll2 = B_coll

all_eqs = eqsJ5 + eqsJ6 + eqsJ7 + eqsJ8 + [eq_coll1, eq_coll2]
vars_list = [a20,a11,a02,b20, a30,a21,a12,a03,b30, a40,a31,a22,a13,a04,b40, a50,a41,a32,a23,a14,a05,b50]

print(f"Total eqs: {len(all_eqs)}, vars: {len(vars_list)}")

# Lambdify
funcs = []
for eq in all_eqs:
    f = sp.lambdify(vars_list, eq, modules=['cmath'])
    funcs.append(f)

def eval_residuals(vals):
    res=[]
    for f in funcs:
        try:
            v = f(*vals)
            res.append(complex(v))
        except:
            res.append(complex(1e6,0))
    return res

def objective(vals):
    res = eval_residuals(vals)
    return sum(abs(r)**2 for r in res)

best_vals = [complex(random.uniform(-1,1), random.uniform(-1,1)) for _ in range(len(vars_list))]
best_obj = objective(best_vals)
print(f"Initial objective: {best_obj}")

iterations=5000
step_size=0.5
for it in range(iterations):
    idx = random.randrange(len(vars_list))
    old_val = best_vals[idx]
    perturb = complex(random.gauss(0, step_size), random.gauss(0, step_size))
    new_vals = best_vals.copy()
    new_vals[idx] = old_val + perturb
    new_obj = objective(new_vals)
    if new_obj < best_obj:
        best_vals = new_vals
        best_obj = new_obj
        if it % 200 ==0:
            print(f"Iter {it}: improved to {best_obj}")
        if best_obj < 1e-12:
            print("Found near-zero!")
            break
    if it % 500 ==0:
        print(f"Iter {it}, best_obj {best_obj}, step {step_size}")
        step_size *= 0.95

print(f"Final best_obj: {best_obj}")
res = eval_residuals(best_vals)
for i,r in enumerate(res):
    if abs(r)>1e-3:
        print(f"eq {i} residual {r}")

if best_obj < 1e-6:
    print("Candidate found!")
else:
    print("No solution found")

print("Done")
