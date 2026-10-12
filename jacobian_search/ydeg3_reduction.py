import sympy as sp, random

x=sp.symbols('x')

# Unknowns
h0,h1, r1_1,r1_2, p1_1,p1_2, p0_2,p0_3, r0_1,r0_2,r0_3, alpha,beta,gamma = sp.symbols('h0 h1 r1_1 r1_2 p1_1 p1_2 p0_2 p0_3 r0_1 r0_2 r0_3 alpha beta gamma')

# Polynomials
h = h0 + h1*x
r1 = 1 + r1_1*x + r1_2*x**2
p1 = p1_1*x + p1_2*x**2
p0 = x + p0_2*x**2 + p0_3*x**3
r0 = r0_1*x + r0_2*x**2 + r0_3*x**3

# Collision constraints:
# p0(1)=0 => 1 + p0_2 + p0_3 =0
# r0(1)=0 => r0_1 + r0_2 + r0_3 =0
eq_coll_p0 = 1 + p0_2 + p0_3
eq_coll_r0 = r0_1 + r0_2 + r0_3

# r2 = alpha * h^2
r2 = alpha * h**2
# p3 = beta * h^3
p3 = beta * h**3
# p2 = (3*beta/(2*alpha)) * h * r1 + gamma * h^2, assume alpha !=0
# To avoid division, multiply: p2 = (3*beta * h * r1 + 2*alpha*gamma * h^2) / (2*alpha)
# Let's define p2 numerator
# p2 = (3*beta*h*r1)/(2*alpha) + gamma*h^2
# For polynomial identity, we need alpha !=0, but we can keep symbolic with division
# Use sympy Rational with alpha in denominator, but for solving we can multiply equations by alpha to clear denominators
# Let's keep p2 as expression with division, sympy will handle
p2 = (3*beta*h*r1)/(2*alpha) + gamma*h**2

# Derivatives
p0p = sp.diff(p0,x)
p1p = sp.diff(p1,x)
r0p = sp.diff(r0,x)
r1p = sp.diff(r1,x)
r2p = sp.diff(r2,x)
p2p = sp.diff(p2,x)
p3p = sp.diff(p3,x)
hprime = sp.diff(h,x)

# Equations y0,y1,y2
# y0: p0' * r1 - p1 * r0' =1
eq_y0 = sp.expand(p0p*r1 - p1*r0p - 1)

# y1: 2 p0' r2 + p1' r1 - p1 r1' -2 p2 r0' =0
eq_y1 = sp.expand(2*p0p*r2 + p1p*r1 - p1*r1p - 2*p2*r0p)

# y2: 2 p1' r2 + p2' r1 - p1 r2' -2 p2 r1' -3 p3 r0' =0
eq_y2 = sp.expand(2*p1p*r2 + p2p*r1 - p1*r2p - 2*p2*r1p - 3*p3*r0p)

print("eq_y0:", eq_y0)
print("eq_y1:", eq_y1)
print("eq_y2:", eq_y2)

# Collect coefficients
poly_y0 = sp.Poly(eq_y0, x)
poly_y1 = sp.Poly(eq_y1, x)
poly_y2 = sp.Poly(eq_y2, x)

# Note eq_y1, eq_y2 have alpha in denominator due to p2, multiply by alpha to clear?
# eq_y1 has term -2*p2*r0p = -2*(3*beta*h*r1/(2*alpha) + gamma*h^2)*r0p = -(3*beta*h*r1/alpha +2*gamma*h^2)*r0p
# So eq_y1 * alpha to clear denominator
eq_y1_cleared = sp.expand(eq_y1 * alpha)
eq_y2_cleared = sp.expand(eq_y2 * alpha)  # p2' also has alpha denominator, p2' = derivative of (3*beta*h*r1/(2*alpha) + gamma*h^2) = (3*beta/(2*alpha))*(h'*r1 + h*r1') + 2*gamma*h*h'
# So eq_y2 * alpha to clear

poly_y1_c = sp.Poly(eq_y1_cleared, x)
poly_y2_c = sp.Poly(eq_y2_cleared, x)

print("\nDegrees:", poly_y0.degree(), poly_y1_c.degree(), poly_y2_c.degree())

eqs = []
# y0 coeffs
for coeff in poly_y0.all_coeffs():
    eqs.append(coeff)
# y1 cleared coeffs
for coeff in poly_y1_c.all_coeffs():
    eqs.append(coeff)
# y2 cleared coeffs
for coeff in poly_y2_c.all_coeffs():
    eqs.append(coeff)
# collision
eqs.append(eq_coll_p0)
eqs.append(eq_coll_r0)

print(f"\nTotal equations: {len(eqs)}")
vars_list = [h0,h1, r1_1,r1_2, p1_1,p1_2, p0_2,p0_3, r0_1,r0_2,r0_3, alpha,beta,gamma]
print(f"Vars: {len(vars_list)}")

# Try to solve with sympy solve? Might be heavy (14 vars, ~17 eqs)
# Let's attempt Groebner or solve for subset

# For numeric search, we need to avoid alpha=0, so set alpha=1 w.l.o.g. by scaling h? Actually r2 = alpha h^2, we can absorb alpha into h scaling: if we scale h by sqrt(alpha), r2 scaling? Let's set alpha=1 to reduce vars, similarly beta can be scaled? Let's attempt to set alpha=1 for numeric search to simplify

print("\n--- Attempt numeric search with alpha=1 ---")
# Set alpha=1
alpha_val=1
# Substitute alpha=1 into equations
eqs_alpha1 = [eq.subs({alpha:1}) for eq in eqs]
vars_alpha1 = [h0,h1, r1_1,r1_2, p1_1,p1_2, p0_2,p0_3, r0_1,r0_2,r0_3, beta,gamma]
print(f"Eqs with alpha=1: {len(eqs_alpha1)}, vars: {len(vars_alpha1)}")

# Try nsolve with random initial guesses
# Need square system: pick 13 equations out of len
# Let's pick first 13
eqs_square = eqs_alpha1[:13]
print(f"Square eqs: {len(eqs_square)}")

for attempt in range(20):
    init = [random.uniform(-2,2) for _ in range(len(vars_alpha1))]
    try:
        sol = sp.nsolve(eqs_square, vars_alpha1, init, tol=1e-14, maxsteps=200, prec=50)
        print(f"\nAttempt {attempt} found solution: {sol}")
        subs = {v:val for v,val in zip(vars_alpha1, sol)}
        # Check remaining eqs
        rem = eqs_alpha1[13:]
        residuals = [complex(eq.subs(subs).evalf()) for eq in rem]
        print("Remaining residuals:", residuals)
        if all(abs(r)<1e-6 for r in residuals):
            print("*** FULL SOLUTION FOUND ***")
            # Build polynomials
            h_val = h.subs(subs)
            r1_val = r1.subs(subs)
            p1_val = p1.subs(subs)
            p0_val = p0.subs(subs)
            r0_val = r0.subs(subs)
            r2_val = r2.subs({alpha:1}).subs(subs)
            p3_val = p3.subs(subs)
            p2_val = p2.subs({alpha:1}).subs(subs)
            print(f"h={h_val}, r1={r1_val}, p1={p1_val}, p0={p0_val}, r0={r0_val}")
            print(f"r2={r2_val}, p3={p3_val}, p2={p2_val}")
            break
    except Exception as e:
        print(f"Attempt {attempt} nsolve failed: {e}")
        continue

print("\n--- Attempt with general alpha, beta, gamma but h linear, r1 quadratic etc. ---")
# For exact solving, try to use Groebner with small degrees? Might be heavy, try with alpha=1, beta=1 to further reduce
# Set alpha=1, beta=1
eqs_alpha1_beta1 = [eq.subs({alpha:1, beta:1}) for eq in eqs]
vars_ab1 = [h0,h1, r1_1,r1_2, p1_1,p1_2, p0_2,p0_3, r0_1,r0_2,r0_3, gamma]
print(f"Eqs alpha=1,beta=1: {len(eqs_alpha1_beta1)}, vars: {len(vars_ab1)}")
eqs_square2 = eqs_alpha1_beta1[:12]
for attempt in range(20):
    init = [random.uniform(-2,2) for _ in range(len(vars_ab1))]
    try:
        sol = sp.nsolve(eqs_square2, vars_ab1, init, tol=1e-14, maxsteps=200, prec=50)
        print(f"\nAttempt {attempt} found solution (alpha=1,beta=1): {sol}")
        subs = {v:val for v,val in zip(vars_ab1, sol)}
        rem = eqs_alpha1_beta1[12:]
        residuals = [complex(eq.subs(subs).evalf()) for eq in rem]
        print("Remaining residuals:", residuals)
        if all(abs(r)<1e-6 for r in residuals):
            print("*** FULL SOLUTION FOUND alpha=1,beta=1 ***")
            break
    except Exception as e:
        print(f"Attempt {attempt} failed: {e}")
        continue

print("\nDone y-degree 3 reduction")
