import sympy as sp, random

x,y=sp.symbols('x y')

# S,T degree 3: 9 coeffs each
v0,v1,w0,w1 = sp.symbols('v0 v1 w0 w1')
s1,s2,s3,s4,s5,s6,s7,s8,s9 = sp.symbols('s1 s2 s3 s4 s5 s6 s7 s8 s9')
t1,t2,t3,t4,t5,t6,t7,t8,t9 = sp.symbols('t1 t2 t3 t4 t5 t6 t7 t8 t9')

A_x0 = -(1+v0)*x**2 + (v0 - v1)*x**3 + v1*x**4
B_x0 = -w0*x**2 + (w0 - w1)*x**3 + w1*x**4

S = s1*x + s2*y + s3*x**2 + s4*x*y + s5*y**2 + s6*x**3 + s7*x**2*y + s8*x*y**2 + s9*y**3
T = t1*x + t2*y + t3*x**2 + t4*x*y + t5*y**2 + t6*x**3 + t7*x**2*y + t8*x*y**2 + t9*y**3

A = A_x0 + y*S
B = B_x0 + y*T

Ax=sp.diff(A,x)
Ay=sp.diff(A,y)
Bx=sp.diff(B,x)
By=sp.diff(B,y)
J=sp.expand(Ax+By+Ax*By-Ay*Bx)
poly=sp.Poly(J, x,y)
eqs=list(poly.as_dict().values())
vars_list = [v0,v1,w0,w1,s1,s2,s3,s4,s5,s6,s7,s8,s9,t1,t2,t3,t4,t5,t6,t7,t8,t9]
print(f"Number eqs: {len(eqs)}, vars: {len(vars_list)}")

funcs = [sp.lambdify(vars_list, eq, modules=['cmath']) for eq in eqs]

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

iterations=10000
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
    print("No solution")

print("Done")
