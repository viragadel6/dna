import sympy as sp, random

x=sp.symbols('x')

h0,h1, r1_1,r1_2, p1_1,p1_2, p0_2,p0_3, r0_1,r0_2,r0_3, beta,gamma = sp.symbols('h0 h1 r1_1 r1_2 p1_1 p1_2 p0_2 p0_3 r0_1 r0_2 r0_3 beta gamma')
# alpha=1
alpha=1

h = h0 + h1*x
r1 = 1 + r1_1*x + r1_2*x**2
p1 = p1_1*x + p1_2*x**2
p0 = x + p0_2*x**2 + p0_3*x**3
r0 = r0_1*x + r0_2*x**2 + r0_3*x**3

r2 = alpha * h**2
p3 = beta * h**3
p2 = (3*beta*h*r1)/(2*alpha) + gamma*h**2

p0p = sp.diff(p0,x)
p1p = sp.diff(p1,x)
r0p = sp.diff(r0,x)
r1p = sp.diff(r1,x)
r2p = sp.diff(r2,x)
p2p = sp.diff(p2,x)
p3p = sp.diff(p3,x)

eq_y0 = sp.expand(p0p*r1 - p1*r0p - 1)
eq_y1 = sp.expand(2*p0p*r2 + p1p*r1 - p1*r1p - 2*p2*r0p)
eq_y2 = sp.expand(2*p1p*r2 + p2p*r1 - p1*r2p - 2*p2*r1p - 3*p3*r0p)

eq_coll_p0 = 1 + p0_2 + p0_3
eq_coll_r0 = r0_1 + r0_2 + r0_3

poly_y0 = sp.Poly(eq_y0, x)
poly_y1 = sp.Poly(eq_y1, x)
poly_y2 = sp.Poly(eq_y2, x)

eqs = []
for c in poly_y0.all_coeffs():
    eqs.append(c)
for c in poly_y1.all_coeffs():
    eqs.append(c)
for c in poly_y2.all_coeffs():
    eqs.append(c)
eqs.append(eq_coll_p0)
eqs.append(eq_coll_r0)

vars_list = [h0,h1, r1_1,r1_2, p1_1,p1_2, p0_2,p0_3, r0_1,r0_2,r0_3, beta,gamma]
print(f"Eqs: {len(eqs)}, vars: {len(vars_list)}")

funcs = [sp.lambdify(vars_list, eq, modules=['cmath']) for eq in eqs]

def eval_residuals(vals):
    res=[]
    for f in funcs:
        try:
            v = f(*vals)
            res.append(complex(v))
        except Exception as e:
            res.append(complex(1e6,0))
    return res

def objective(vals):
    res = eval_residuals(vals)
    return sum(abs(r)**2 for r in res)

best_vals = [complex(random.uniform(-2,2), random.uniform(-2,2)) for _ in range(len(vars_list))]
best_obj = objective(best_vals)
print(f"Initial obj: {best_obj}")

iterations=20000
step_size=0.5
for it in range(iterations):
    idx = random.randrange(len(vars_list))
    old = best_vals[idx]
    perturb = complex(random.gauss(0, step_size), random.gauss(0, step_size))
    new_vals = best_vals.copy()
    new_vals[idx] = old + perturb
    new_obj = objective(new_vals)
    if new_obj < best_obj:
        best_vals = new_vals
        best_obj = new_obj
        if it % 500 ==0:
            print(f"Iter {it}: improved to {best_obj}")
        if best_obj < 1e-12:
            print("Found near-zero!")
            break
    if it % 1000 ==0:
        print(f"Iter {it}, best_obj {best_obj}, step {step_size}")
        step_size *= 0.95

print(f"Final best_obj: {best_obj}")
res = eval_residuals(best_vals)
for i,r in enumerate(res):
    if abs(r)>1e-3:
        print(f"eq {i} residual {r}")

if best_obj < 1e-6:
    print("Candidate found!")
    subs = {v:val for v,val in zip(vars_list, best_vals)}
    print(f"h={h.subs(subs)}, r1={r1.subs(subs)}, p1={p1.subs(subs)}, p0={p0.subs(subs)}, r0={r0.subs(subs)}")
    print(f"r2={r2.subs(subs)}, p3={p3.subs(subs)}, p2={p2.subs(subs)}")
    # Build P,R
    # P = p0 + p1 y + p2 y^2 + p3 y^3
    # R = r0 + r1 y + r2 y^2
    # Then Q = R + lambda P, with lambda = -r0_1 (since r0_1 = -lambda)
    # For our case lambda = -r0_1
    # Let's compute lambda
    lam = -r0_1
    lam_val = lam.subs(subs)
    print(f"lambda={lam_val}")
    # Build full F1,F2
    y=sp.symbols('y')
    P_full = p0.subs(subs) + p1.subs(subs)*y + p2.subs(subs)*y**2 + p3.subs(subs)*y**3
    R_full = r0.subs(subs) + r1.subs(subs)*y + r2.subs(subs)*y**2
    Q_full = R_full + lam_val * P_full
    print(f"P={sp.expand(P_full)}")
    print(f"Q={sp.expand(Q_full)}")
    # Check collision: P(0,0)=0, P(1,0)=0?
    print(f"P(0,0)={P_full.subs({x:0})}, P(1,0) with y=0: {p0.subs(subs).subs({x:1})}")
    print(f"Q(0,0)={R_full.subs({x:0})}, Q(1,0)={r0.subs(subs).subs({x:1}) + lam_val * p0.subs(subs).subs({x:1})}")
else:
    print("No solution found")

print("Done")
