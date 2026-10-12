import sympy as sp, random, math, cmath

x,y=sp.symbols('x y')

# Setup d=4 full system as before but with 15 vars
a20,a11,a02,b20, a30,a21,a12,a03,b30, a40,a31,a22,a13,a04,b40 = sp.symbols('a20 a11 a02 b20 a30 a21 a12 a03 b30 a40 a31 a22 a13 a04 b40')

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

A4 = a40*x**4 + a31*x**3*y + a22*x**2*y**2 + a13*x*y**3 + a04*y**4
# B4 generic with b31,b22,b13,b04 expressed from J3
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

JacA2B4 = sp.diff(A2,x)*sp.diff(B4,y) - sp.diff(A2,y)*sp.diff(B4,x)
JacA3B3 = sp.diff(A3,x)*sp.diff(B3,y) - sp.diff(A3,y)*sp.diff(B3,x)
JacA4B2 = sp.diff(A4,x)*sp.diff(B2,y) - sp.diff(A4,y)*sp.diff(B2,x)
J4 = sp.expand(JacA2B4 + JacA3B3 + JacA4B2)
polyJ4 = sp.Poly(J4, x,y)
eqsJ4 = list(polyJ4.as_dict().values())

JacA3B4 = sp.diff(A3,x)*sp.diff(B4,y) - sp.diff(A3,y)*sp.diff(B4,x)
JacA4B3 = sp.diff(A4,x)*sp.diff(B3,y) - sp.diff(A4,y)*sp.diff(B3,x)
J5 = sp.expand(JacA3B4 + JacA4B3)
polyJ5 = sp.Poly(J5, x,y)
eqsJ5 = list(polyJ5.as_dict().values())

JacA4B4 = sp.diff(A4,x)*sp.diff(B4,y) - sp.diff(A4,y)*sp.diff(B4,x)
J6 = sp.expand(JacA4B4)
polyJ6 = sp.Poly(J6, x,y)
eqsJ6 = list(polyJ6.as_dict().values())

A_coll = (A2 + A3 + A4).subs({x:1,y:0})
B_coll = (B2 + B3 + B4).subs({x:1,y:0})
eq_coll1 = A_coll + 1
eq_coll2 = B_coll

all_eqs = eqsJ4 + eqsJ5 + eqsJ6 + [eq_coll1, eq_coll2]
print(f"Total eqs: {len(all_eqs)}")
vars_list = [a20,a11,a02,b20, a30,a21,a12,a03,b30, a40,a31,a22,a13,a04,b40]
print(f"Vars: {len(vars_list)}")

# We'll attempt random search with hill climbing
# Objective: sum |eq|^2

def eval_residuals(vals):
    # vals: list of complex numbers
    # Evaluate each eq
    res=[]
    # Use subs for sympy eval (slow but okay for small)
    subs_dict = {var: val for var,val in zip(vars_list, vals)}
    for eq in all_eqs:
        try:
            v = complex(eq.subs(subs_dict).evalf())
        except Exception as e:
            v = complex(1e6,0)
        res.append(v)
    return res

def objective(vals):
    res = eval_residuals(vals)
    return sum(abs(r)**2 for r in res)

# Random hill climbing
best_vals = [complex(random.uniform(-1,1), random.uniform(-1,1)) for _ in range(len(vars_list))]
best_obj = objective(best_vals)
print(f"Initial objective: {best_obj}")

# Simple simulated annealing
# For iteration in range
iterations=2000
step_size=0.5
for it in range(iterations):
    # Perturb one variable
    idx = random.randrange(len(vars_list))
    old_val = best_vals[idx]
    # Perturb with complex gaussian
    perturb = complex(random.gauss(0, step_size), random.gauss(0, step_size))
    new_vals = best_vals.copy()
    new_vals[idx] = old_val + perturb
    new_obj = objective(new_vals)
    if new_obj < best_obj:
        best_vals = new_vals
        best_obj = new_obj
        print(f"Iter {it}: improved objective to {best_obj}")
        if best_obj < 1e-12:
            print("Found near-zero solution!")
            break
    # Occasionally reduce step size
    if it % 200 ==0:
        step_size *= 0.9
        print(f"Iter {it}, best_obj {best_obj}, step_size {step_size}")

print(f"Final best objective: {best_obj}")
print(f"Best vals: {best_vals}")
res = eval_residuals(best_vals)
for eq, r in zip(all_eqs, res):
    print(f"eq residual: {r}")

# Check if near zero, then we have candidate
if best_obj < 1e-6:
    # Build A,B
    subs_dict = {var: val for var,val in zip(vars_list, best_vals)}
    A_full = (A2 + A3 + A4).subs(subs_dict)
    B_full = (B2 + B3 + B4).subs(subs_dict)
    print("A_full:", sp.expand(A_full))
    print("B_full:", sp.expand(B_full))
    # Check collision
    print("A(1,0):", A_full.subs({x:1,y:0}))
    print("B(1,0):", B_full.subs({x:1,y:0}))
    # Check Jacobian
    Ax = sp.diff(A_full,x)
    Ay = sp.diff(A_full,y)
    Bx = sp.diff(B_full,x)
    By = sp.diff(B_full,y)
    J = sp.expand(Ax+By+Ax*By-Ay*Bx)
    print("J:", J)

print("Done random search")
