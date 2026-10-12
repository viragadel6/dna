import sympy as sp
import argparse
import heapq
import itertools
import os
import sys
import traceback


def emit(*values):
    print(*values, flush=True)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def nonnegative_integer(text):
    try:
        value = int(text)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("The degree must be an integer.") from exc
    if value < 0:
        raise argparse.ArgumentTypeError("The degree must be nonnegative.")
    return value


def parse_degree(argv=None):
    parser = argparse.ArgumentParser(
        description="Exact bounded-degree search for a normalized two-dimensional Jacobian counterexample.",
        epilog=(
            "Large degree bounds can exceed available time and memory. "
            "Excluding one degree bound does not prove the Jacobian conjecture."
        ),
    )
    parser.add_argument(
        "degree",
        nargs="?",
        type=nonnegative_integer,
        help="Maximum total degree of the nonlinear parts P and Q.",
    )
    parser.add_argument(
        "-d",
        "--degree",
        dest="option_degree",
        type=nonnegative_integer,
        help="Maximum total degree of the nonlinear parts P and Q.",
    )
    arguments = parser.parse_args(argv)
    if arguments.degree is not None and arguments.option_degree is not None:
        parser.error("Specify the degree only once.")
    degree = (
        arguments.degree
        if arguments.degree is not None
        else arguments.option_degree
    )
    if degree is None:
        try:
            text = input("Maximum nonlinear degree d: ")
        except EOFError:
            parser.error("No degree was supplied.")
        try:
            degree = nonnegative_integer(text)
        except argparse.ArgumentTypeError as exc:
            parser.error(str(exc))
    return degree


def all_monomials(deg_min, deg_max):
    monomials = []
    for total_degree in range(deg_min, deg_max + 1):
        for i in range(total_degree + 1):
            j = total_degree - i
            monomials.append((i, j))
    return tuple(monomials)


def monomial_expression(generators, exponents):
    result = sp.S.One
    for generator, exponent in zip(generators, exponents):
        result *= generator ** exponent
    return result


def build_poly_with_vars(name, x, y, monomials):
    coefficients = []
    terms = []
    for i, j in monomials:
        coefficient = sp.Symbol(f"{name}_{i}_{j}")
        coefficients.append(coefficient)
        terms.append(coefficient * x ** i * y ** j)
    return sp.expand(sp.Add(*terms)), tuple(coefficients)


def rational_groebner(equations, generators):
    generators = tuple(generators)
    require(bool(generators), "A polynomial Groebner computation needs generators.")
    return sp.groebner(
        [sp.expand(equation) for equation in equations],
        *generators,
        order="lex",
        domain=sp.QQ,
        method="buchberger",
    )


def is_unit_basis(basis):
    return len(basis.polys) == 1 and basis.polys[0].as_expr() == sp.S.One


def print_basis(title, basis):
    emit(title)
    emit("Coefficient domain:", basis.domain)
    emit("Lexicographic generator order:", basis.gens)
    emit([sp.expand(polynomial.as_expr()) for polynomial in basis.polys])


def leading_exponent(polynomial):
    require(not polynomial.is_zero, "A zero polynomial has no leading monomial.")
    return tuple(polynomial.monoms(order="lex")[0])


def exponent_divides(divisor, dividend):
    return all(a <= b for a, b in zip(divisor, dividend))


def divide_with_polynomials(polynomial, divisor_expressions, generators):
    if not divisor_expressions:
        return [], polynomial
    quotients, remainder = sp.reduced(
        polynomial.as_expr(),
        divisor_expressions,
        *generators,
        order="lex",
        domain=sp.QQ,
    )
    quotient_polynomials = [
        sp.Poly(quotient, *generators, domain=sp.QQ)
        for quotient in quotients
    ]
    remainder_polynomial = sp.Poly(
        remainder, *generators, domain=sp.QQ
    )
    return quotient_polynomials, remainder_polynomial


def add_certificate_multiple(destination, source, multiplier, zero):
    if multiplier.is_zero:
        return
    for index, coefficient in source.items():
        value = destination.get(index, zero) + multiplier * coefficient
        if value.is_zero:
            destination.pop(index, None)
        else:
            destination[index] = value


def recover_membership_certificates(equations, target_basis):
    generators = tuple(target_basis.gens)
    originals = [
        sp.Poly(equation, *generators, domain=sp.QQ)
        for equation in equations
    ]
    targets = list(target_basis.polys)
    if not targets:
        return []

    zero = sp.Poly(0, *generators, domain=sp.QQ)
    one = sp.Poly(1, *generators, domain=sp.QQ)
    basis_polynomials = []
    basis_expressions = []
    basis_certificates = []
    leading_exponents = []
    critical_pairs = []

    def append_generator(polynomial, certificate):
        require(not polynomial.is_zero, "Cannot append a zero basis generator.")
        scale = sp.S.One / polynomial.LC()
        normalized = polynomial.monic()
        normalized_certificate = {
            index: coefficient * scale
            for index, coefficient in certificate.items()
            if not coefficient.is_zero
        }
        new_index = len(basis_polynomials)
        new_leading = leading_exponent(normalized)
        for old_index, old_leading in enumerate(leading_exponents):
            common = tuple(
                max(a, b) for a, b in zip(old_leading, new_leading)
            )
            heapq.heappush(
                critical_pairs,
                (sum(common), common, old_index, new_index),
            )
        basis_polynomials.append(normalized)
        basis_expressions.append(normalized.as_expr())
        basis_certificates.append(normalized_certificate)
        leading_exponents.append(new_leading)

    def target_certificates_if_available():
        if not basis_expressions:
            return None
        target_quotients = []
        for target in targets:
            quotients, remainder = divide_with_polynomials(
                target, basis_expressions, generators
            )
            if not remainder.is_zero:
                return None
            target_quotients.append(quotients)

        certificates = []
        for quotients in target_quotients:
            certificate = {}
            for quotient, source in zip(quotients, basis_certificates):
                add_certificate_multiple(
                    certificate, source, quotient, zero
                )
            certificates.append(certificate)
        return certificates

    for original_index, original in enumerate(originals):
        if original.is_zero:
            continue
        quotients, remainder = divide_with_polynomials(
            original, basis_expressions, generators
        )
        if remainder.is_zero:
            continue
        certificate = {original_index: one}
        for quotient, source in zip(quotients, basis_certificates):
            add_certificate_multiple(
                certificate, source, -quotient, zero
            )
        append_generator(remainder, certificate)
        if basis_polynomials[-1].as_expr() == sp.S.One:
            break

    certificates = target_certificates_if_available()
    if certificates is not None:
        return certificates

    while critical_pairs:
        _, common, left_index, right_index = heapq.heappop(critical_pairs)
        left_difference = tuple(
            a - b
            for a, b in zip(common, leading_exponents[left_index])
        )
        right_difference = tuple(
            a - b
            for a, b in zip(common, leading_exponents[right_index])
        )
        left_multiplier = sp.Poly(
            monomial_expression(generators, left_difference),
            *generators,
            domain=sp.QQ,
        )
        right_multiplier = sp.Poly(
            monomial_expression(generators, right_difference),
            *generators,
            domain=sp.QQ,
        )
        s_polynomial = (
            left_multiplier * basis_polynomials[left_index]
            - right_multiplier * basis_polynomials[right_index]
        )
        quotients, remainder = divide_with_polynomials(
            s_polynomial, basis_expressions, generators
        )
        if remainder.is_zero:
            continue

        certificate = {}
        add_certificate_multiple(
            certificate,
            basis_certificates[left_index],
            left_multiplier,
            zero,
        )
        add_certificate_multiple(
            certificate,
            basis_certificates[right_index],
            -right_multiplier,
            zero,
        )
        for quotient, source in zip(quotients, basis_certificates):
            add_certificate_multiple(
                certificate, source, -quotient, zero
            )
        append_generator(remainder, certificate)

        certificates = target_certificates_if_available()
        if certificates is not None:
            return certificates

    raise ArithmeticError(
        "Buchberger membership reconstruction finished without certifying "
        "the claimed Groebner basis."
    )


def certified_basis_division(expression, basis, label):
    expression = sp.expand(expression)
    if basis.polys:
        quotients, remainder = basis.reduce(expression)
    else:
        quotients, remainder = [], expression
    quotients = [sp.expand(quotient) for quotient in quotients]
    remainder = sp.expand(remainder)
    reconstructed = remainder
    for index, quotient in enumerate(quotients):
        reconstructed += quotient * basis.polys[index].as_expr()
        emit(f"{label}: quotient for G_{index + 1} =", quotient)
    identity_residual = sp.expand(expression - reconstructed)
    require(
        identity_residual == 0,
        f"The exact division identity failed for {label}.",
    )
    emit(f"{label}: remainder =", remainder)
    emit(f"{label}: division identity residual =", identity_residual)
    return remainder


def certify_reduced_groebner(basis, equations):
    require(basis.domain == sp.QQ, "The basis is not over QQ.")
    require(str(basis.order) == "lex", "The basis order is not lexicographic.")
    polynomials = list(basis.polys)
    leading = []

    for index, polynomial in enumerate(polynomials):
        require(not polynomial.is_zero, "The basis contains zero.")
        require(
            polynomial.LC() == sp.S.One,
            f"G_{index + 1} is not monic.",
        )
        leading.append(leading_exponent(polynomial))

    for index, polynomial in enumerate(polynomials):
        for exponent, coefficient in polynomial.terms(order="lex"):
            if coefficient == 0:
                continue
            for divisor_index, divisor in enumerate(leading):
                if index == divisor_index and tuple(exponent) == leading[index]:
                    continue
                require(
                    not exponent_divides(divisor, exponent),
                    f"G_{index + 1} is not reduced with respect to the basis.",
                )

    emit("Monicity and reducedness checks:", True)

    for index, equation in enumerate(equations, start=1):
        remainder = certified_basis_division(
            equation, basis, f"E_{index}"
        )
        require(
            remainder == 0,
            f"E_{index} is not in the ideal generated by the reported basis.",
        )

    for left_index in range(len(polynomials)):
        for right_index in range(left_index + 1, len(polynomials)):
            common = tuple(
                max(a, b)
                for a, b in zip(leading[left_index], leading[right_index])
            )
            left_difference = tuple(
                a - b for a, b in zip(common, leading[left_index])
            )
            right_difference = tuple(
                a - b for a, b in zip(common, leading[right_index])
            )
            s_polynomial = sp.expand(
                monomial_expression(basis.gens, left_difference)
                * polynomials[left_index].as_expr()
                - monomial_expression(basis.gens, right_difference)
                * polynomials[right_index].as_expr()
            )
            label = f"S_{left_index + 1}_{right_index + 1}"
            emit(f"{label} =", s_polynomial)
            remainder = certified_basis_division(
                s_polynomial, basis, label
            )
            require(
                remainder == 0,
                f"Buchberger's criterion failed for {label}.",
            )

    emit("All Buchberger critical-pair checks:", True)
    emit("Reconstructing exact membership certificates in the original equations.")
    certificates = recover_membership_certificates(equations, basis)
    generators = tuple(basis.gens)
    zero = sp.Poly(0, *generators, domain=sp.QQ)
    originals = [
        sp.Poly(equation, *generators, domain=sp.QQ)
        for equation in equations
    ]

    for basis_index, polynomial in enumerate(polynomials):
        certificate = certificates[basis_index]
        reconstructed = zero
        for equation_index, original in enumerate(originals):
            multiplier = certificate.get(equation_index, zero)
            emit(
                f"C_{basis_index + 1}_{equation_index + 1} =",
                multiplier.as_expr(),
            )
            reconstructed += multiplier * original
        residual = polynomial - reconstructed
        require(
            residual.is_zero,
            f"The original-ideal membership certificate failed for "
            f"G_{basis_index + 1}.",
        )
        emit(
            f"G_{basis_index + 1} - "
            f"sum(C_{basis_index + 1}_j * E_j) =",
            residual.as_expr(),
        )

    emit("Exact reduced Groebner basis and equality of ideals certified:", True)


def certify_constant_system(equations):
    values = [sp.Rational(equation) for equation in equations]
    nonzero_index = next(
        (index for index, value in enumerate(values) if value != 0),
        None,
    )
    require(
        nonzero_index is not None,
        "The coefficient-free normalized collision system was unexpectedly consistent.",
    )
    emit("Reduced Groebner basis over QQ:")
    emit([sp.S.One])
    reconstructed = sp.S.Zero
    for index, value in enumerate(values):
        multiplier = (
            sp.S.One / values[nonzero_index]
            if index == nonzero_index
            else sp.S.Zero
        )
        emit(f"C_1_{index + 1} =", multiplier)
        reconstructed += multiplier * value
        emit(f"E_{index + 1}: quotient for G_1 =", value)
        emit(f"E_{index + 1}: remainder =", sp.S.Zero)
    residual = sp.expand(sp.S.One - reconstructed)
    require(residual == 0, "The constant unit-ideal certificate failed.")
    emit("1 - sum(C_1_j * E_j) =", residual)
    emit("Exact unit-ideal certificate verified:", True)


def missing_pure_power_indices(basis):
    present = set()
    for polynomial in basis.polys:
        exponent = leading_exponent(polynomial)
        support = [
            index for index, power in enumerate(exponent) if power != 0
        ]
        if len(support) == 1:
            present.add(support[0])
    return [
        index for index in range(len(basis.gens))
        if index not in present
    ]


def integer_sequence():
    yield sp.S.Zero
    for magnitude in itertools.count(1):
        value = sp.Integer(magnitude)
        yield value
        yield -value


def obtain_zero_dimensional_slice(basis):
    require(not is_unit_basis(basis), "Cannot sample a unit ideal.")
    generators = tuple(basis.gens)
    current = basis
    slice_count = 0

    while not current.is_zero_dimensional:
        missing = missing_pure_power_indices(current)
        require(
            bool(missing),
            "The dimension test and leading-monomial test disagree.",
        )
        variable = generators[missing[-1]]
        old_missing_count = len(missing)
        current_equations = [
            polynomial.as_expr() for polynomial in current.polys
        ]

        for value in integer_sequence():
            slice_equation = sp.expand(variable - value)
            trial = rational_groebner(
                current_equations + [slice_equation],
                generators,
            )
            emit("Rational slice trial:", variable, "=", value)
            print_basis("Rational-slice trial basis:", trial)
            emit("Trial ideal is the unit ideal:", is_unit_basis(trial))
            if is_unit_basis(trial):
                continue

            new_missing_count = len(missing_pure_power_indices(trial))
            require(
                new_missing_count < old_missing_count,
                "A successful rational slice did not advance the dimension reduction.",
            )
            require(
                sp.expand(trial.reduce(slice_equation)[1]) == 0,
                "The selected slice equation is absent from the new ideal.",
            )
            slice_count += 1
            require(
                slice_count <= len(generators),
                "Too many coordinate slices were required.",
            )
            emit(f"Accepted slice {slice_count}:", slice_equation, "= 0")
            current = trial
            break

    require(
        not is_unit_basis(current) and current.is_zero_dimensional,
        "A proper zero-dimensional slice was not obtained.",
    )
    print_basis("Proper zero-dimensional slice basis:", current)
    return current


def univariate_eliminant(basis, variable):
    require(
        basis.gens[-1] == variable and str(basis.order) == "lex",
        "The requested variable is not last in the lexicographic elimination order.",
    )
    candidates = []
    for polynomial in basis.polys:
        expression = polynomial.as_expr()
        if expression.free_symbols.issubset({variable}):
            candidate = sp.Poly(expression, variable, domain=sp.QQ)
            if not candidate.is_zero:
                require(
                    candidate.degree() > 0,
                    "A nonzero constant appeared in a proper elimination ideal.",
                )
                candidates.append(candidate.monic())
    require(
        bool(candidates),
        "A zero-dimensional ideal did not provide a univariate eliminant.",
    )
    eliminant = min(candidates, key=lambda polynomial: polynomial.degree())
    for candidate in candidates:
        require(
            candidate.rem(eliminant).is_zero,
            "The univariate elimination polynomials are inconsistent.",
        )
    return eliminant


def radicalize_zero_dimensional_ideal(basis):
    require(
        not is_unit_basis(basis) and basis.is_zero_dimensional,
        "Radicalization requires a proper zero-dimensional ideal.",
    )
    generators = tuple(basis.gens)
    equations = [polynomial.as_expr() for polynomial in basis.polys]
    squarefree_constraints = []

    for index, variable in enumerate(generators, start=1):
        reordered = tuple(
            generator for generator in generators if generator != variable
        ) + (variable,)
        elimination_basis = rational_groebner(equations, reordered)
        require(
            not is_unit_basis(elimination_basis)
            and elimination_basis.is_zero_dimensional,
            "A coordinate elimination unexpectedly lost the proper finite variety.",
        )
        print_basis(
            f"Coordinate-elimination basis {index}:",
            elimination_basis,
        )
        eliminant = univariate_eliminant(elimination_basis, variable)
        squarefree = eliminant.sqf_part().monic()
        require(squarefree.degree() > 0, "A squarefree constraint became constant.")
        require(
            eliminant.rem(squarefree).is_zero,
            "The squarefree part does not divide the eliminant.",
        )
        require(
            squarefree.gcd(squarefree.diff()).degree() == 0,
            "A purported squarefree polynomial is not squarefree.",
        )
        emit("Elimination variable:", variable)
        emit("Coordinate eliminant:", eliminant.as_expr())
        emit("Squarefree coordinate constraint:", squarefree.as_expr(), "= 0")
        squarefree_constraints.append(squarefree.as_expr())
        equations = [
            polynomial.as_expr()
            for polynomial in elimination_basis.polys
        ] + [squarefree.as_expr()]

    radical_basis = rational_groebner(equations, generators)
    require(
        not is_unit_basis(radical_basis)
        and radical_basis.is_zero_dimensional,
        "The radical construction did not retain a proper finite variety.",
    )
    for constraint in squarefree_constraints:
        require(
            sp.expand(radical_basis.reduce(constraint)[1]) == 0,
            "A squarefree coordinate constraint was lost.",
        )
    print_basis("Reduced basis of the zero-dimensional radical ideal:", radical_basis)
    return radical_basis


def obtain_shape_representation(radical_basis):
    generators = tuple(radical_basis.gens)
    require(bool(generators), "A shape representation needs coefficient variables.")
    require(
        not is_unit_basis(radical_basis)
        and radical_basis.is_zero_dimensional,
        "A shape representation requires a proper finite variety.",
    )
    theta = sp.Symbol("theta")
    require(
        theta not in generators,
        "The primitive-element symbol collides with an existing generator.",
    )
    radical_equations = [
        polynomial.as_expr() for polynomial in radical_basis.polys
    ]

    for integer_parameter in itertools.count(0):
        parameter = sp.Integer(integer_parameter)
        linear_form = sp.expand(
            sp.Add(
                *[
                    parameter ** index * variable
                    for index, variable in enumerate(generators)
                ]
            )
        )
        primitive_equation = sp.expand(theta - linear_form)
        shape_basis = rational_groebner(
            radical_equations + [primitive_equation],
            generators + (theta,),
        )
        require(
            not is_unit_basis(shape_basis)
            and shape_basis.is_zero_dimensional,
            "Adjoining a linear primitive element changed the finite variety incorrectly.",
        )
        emit("Primitive-element trial parameter:", parameter)
        emit("Primitive-element linear form:", linear_form)
        print_basis("Shape-position trial basis:", shape_basis)
        eliminant = univariate_eliminant(shape_basis, theta)

        representatives = {}
        for variable in generators:
            _, remainder = shape_basis.reduce(variable)
            representatives[variable] = sp.expand(remainder)

        shape_condition = all(
            representative.free_symbols.issubset({theta})
            for representative in representatives.values()
        )
        emit("Shape criterion satisfied:", shape_condition)
        if not shape_condition:
            continue

        require(
            eliminant.gcd(eliminant.diff()).degree() == 0,
            "The radical shape eliminant is not squarefree.",
        )
        representatives = {
            variable: sp.Poly(
                representative, theta, domain=sp.QQ
            ).rem(eliminant).as_expr()
            for variable, representative in representatives.items()
        }

        for index, equation in enumerate(
            radical_equations + [primitive_equation],
            start=1,
        ):
            substituted = sp.expand(equation.xreplace(representatives))
            remainder = sp.Poly(
                substituted, theta, domain=sp.QQ
            ).rem(eliminant)
            emit(f"Shape identity {index}: remainder =", remainder.as_expr())
            require(
                remainder.is_zero,
                "The shape representation failed an exact defining equation.",
            )

        emit("Shape eliminant:", eliminant.as_expr())
        for variable in generators:
            emit(
                f"{variable} modulo the shape eliminant =",
                representatives[variable],
            )
        return theta, eliminant, representatives

    raise ArithmeticError("Primitive-element enumeration terminated unexpectedly.")


def construct_exact_algebraic_point(basis):
    finite_basis = obtain_zero_dimensional_slice(basis)
    radical_basis = radicalize_zero_dimensional_ideal(finite_basis)
    theta, eliminant, representatives = obtain_shape_representation(radical_basis)

    content, factors = eliminant.factor_list()
    require(bool(factors), "The shape eliminant has no nonconstant factor.")
    factorization_expression = sp.sympify(content)
    emit("Shape-eliminant factorization content:", content)
    for index, (factor, multiplicity) in enumerate(factors, start=1):
        emit(f"Factor {index}:", factor.as_expr())
        emit(f"Factor {index} multiplicity:", multiplicity)
        factorization_expression *= factor.as_expr() ** multiplicity
    require(
        sp.Poly(
            sp.expand(factorization_expression),
            theta,
            domain=sp.QQ,
        ) == eliminant,
        "The exact univariate factorization identity failed.",
    )

    defining_polynomial = factors[0][0].monic()
    require(
        defining_polynomial.degree() > 0,
        "The selected root-defining polynomial is constant.",
    )
    require(
        eliminant.rem(defining_polynomial).is_zero,
        "The selected defining polynomial does not divide the eliminant.",
    )

    if defining_polynomial.degree() == 1:
        alpha = (
            -defining_polynomial.nth(0) / defining_polynomial.nth(1)
        )
    else:
        alpha = sp.CRootOf(
            defining_polynomial.as_expr(),
            theta,
            0,
        )

    require(
        not alpha.free_symbols,
        "The selected algebraic number has free parameters.",
    )
    computed_minimal_polynomial = sp.Poly(
        sp.minimal_polynomial(alpha, theta),
        theta,
        domain=sp.QQ,
    ).monic()
    require(
        computed_minimal_polynomial == defining_polynomial,
        "The selected algebraic number has the wrong minimal polynomial.",
    )

    reduced_representatives = {
        variable: sp.Poly(
            representative, theta, domain=sp.QQ
        ).rem(defining_polynomial).as_expr()
        for variable, representative in representatives.items()
    }

    emit("Selected defining polynomial f(theta):", defining_polynomial.as_expr())
    emit("Selected root index in SymPy CRootOf ordering:", 0)
    emit("alpha =", alpha)
    emit("Computed monic minimal polynomial of alpha:", computed_minimal_polynomial.as_expr())
    emit("Exact arithmetic for the verification: QQ[theta]/(f(theta)), theta = alpha.")
    return theta, defining_polynomial, alpha, reduced_representatives


def divide_univariate_exact(expression, theta, defining_polynomial):
    polynomial = sp.Poly(sp.expand(expression), theta, domain=sp.QQ)
    quotient, remainder = polynomial.div(defining_polynomial)
    require(
        (
            polynomial
            - defining_polynomial * quotient
            - remainder
        ).is_zero,
        "An exact univariate division identity failed.",
    )
    return quotient.as_expr(), remainder.as_expr()


def divide_bivariate_coefficients_exact(
    expression, x, y, theta, defining_polynomial
):
    expression = sp.expand(expression)
    coefficient_ring = sp.QQ.poly_ring(theta)
    polynomial = sp.Poly(
        expression, x, y, domain=coefficient_ring
    )
    quotient_terms = []
    remainder_terms = []
    for (i, j), coefficient in polynomial.terms():
        quotient, remainder = divide_univariate_exact(
            coefficient, theta, defining_polynomial
        )
        spatial_monomial = x ** i * y ** j
        quotient_terms.append(quotient * spatial_monomial)
        remainder_terms.append(remainder * spatial_monomial)
    quotient_expression = sp.expand(sp.Add(*quotient_terms))
    remainder_expression = sp.expand(sp.Add(*remainder_terms))
    residual = sp.expand(
        expression
        - defining_polynomial.as_expr() * quotient_expression
        - remainder_expression
    )
    require(
        residual == 0,
        "Coefficientwise algebraic reduction failed its polynomial identity.",
    )
    return quotient_expression, remainder_expression


def evaluate_bivariate_polynomial(polynomial, point):
    terms = []
    for (i, j), coefficient in polynomial.terms():
        terms.append(coefficient * point[0] ** i * point[1] ** j)
    return sp.expand(sp.Add(*terms))


def verify_and_print_solution(
    degree,
    x,
    y,
    F1,
    F2,
    unknowns,
    equations,
    theta,
    defining_polynomial,
    alpha,
    representatives,
):
    require(
        set(representatives) == set(unknowns),
        "The coefficient assignment is incomplete.",
    )
    for variable in unknowns:
        representative = sp.Poly(
            representatives[variable], theta, domain=sp.QQ
        )
        require(
            representative.degree() < defining_polynomial.degree(),
            f"The representative for {variable} is not reduced.",
        )

    root_replacement = {theta: alpha}
    actual_coefficients = {
        variable: sp.expand(
            representatives[variable].xreplace(root_replacement)
        )
        for variable in unknowns
    }

    emit("Explicit coefficient assignment:")
    for variable in unknowns:
        emit(f"{variable} =", actual_coefficients[variable])

    F1_representative = sp.expand(F1.xreplace(representatives))
    F2_representative = sp.expand(F2.xreplace(representatives))
    require(
        F1_representative.free_symbols.issubset({x, y, theta})
        and F2_representative.free_symbols.issubset({x, y, theta}),
        "Unassigned coefficient symbols remain in the polynomial map.",
    )

    P_representative = sp.expand(F1_representative - x)
    Q_representative = sp.expand(F2_representative - y)
    F1_actual = sp.expand(F1_representative.xreplace(root_replacement))
    F2_actual = sp.expand(F2_representative.xreplace(root_replacement))
    P_actual = sp.expand(P_representative.xreplace(root_replacement))
    Q_actual = sp.expand(Q_representative.xreplace(root_replacement))

    require(
        F1_actual.free_symbols.issubset({x, y})
        and F2_actual.free_symbols.issubset({x, y}),
        "The explicit map has nonconstant, unresolved coefficient parameters.",
    )

    emit("Explicit P(x,y):", P_actual)
    emit("Explicit Q(x,y):", Q_actual)
    emit("Explicit F1(x,y):", F1_actual)
    emit("Explicit F2(x,y):", F2_actual)
    emit("F1 representative in QQ[theta,x,y]:", F1_representative)
    emit("F2 representative in QQ[theta,x,y]:", F2_representative)

    for index, equation in enumerate(equations, start=1):
        substituted = sp.expand(equation.xreplace(representatives))
        quotient, remainder = divide_univariate_exact(
            substituted, theta, defining_polynomial
        )
        emit(f"E_{index} after coefficient substitution:", substituted)
        emit(f"E_{index} quotient by f(theta):", quotient)
        emit(f"E_{index} remainder by f(theta):", remainder)
        require(
            remainder == 0,
            f"The coefficient assignment violates E_{index}.",
        )

    coefficient_ring = sp.QQ.poly_ring(theta)
    for label, nonlinear_part in (
        ("P", P_representative),
        ("Q", Q_representative),
    ):
        polynomial = sp.Poly(
            nonlinear_part, x, y, domain=coefficient_ring
        )
        support = [
            exponent
            for exponent, coefficient in polynomial.terms()
            if coefficient != 0
        ]
        support_valid = all(
            2 <= sum(exponent) <= degree for exponent in support
        )
        emit(f"{label} monomial support:", support)
        emit(f"{label} belongs to (x,y)^2 and has degree at most d:", support_valid)
        require(
            support_valid,
            f"The nonlinear part {label} violates the normalized degree support.",
        )

    representative_dF1x = sp.diff(F1_representative, x)
    representative_dF1y = sp.diff(F1_representative, y)
    representative_dF2x = sp.diff(F2_representative, x)
    representative_dF2y = sp.diff(F2_representative, y)
    determinant_representative = sp.expand(
        representative_dF1x * representative_dF2y
        - representative_dF1y * representative_dF2x
        - 1
    )

    actual_dF1x = sp.diff(F1_actual, x)
    actual_dF1y = sp.diff(F1_actual, y)
    actual_dF2x = sp.diff(F2_actual, x)
    actual_dF2y = sp.diff(F2_actual, y)
    actual_determinant_raw = sp.expand(
        actual_dF1x * actual_dF2y
        - actual_dF1y * actual_dF2x
        - 1
    )
    require(
        sp.expand(
            actual_determinant_raw
            - determinant_representative.xreplace(root_replacement)
        ) == 0,
        "Differentiation of the explicit map disagrees with its exact representation.",
    )

    determinant_quotient, determinant_remainder = (
        divide_bivariate_coefficients_exact(
            determinant_representative,
            x,
            y,
            theta,
            defining_polynomial,
        )
    )
    emit("Expanded det(JF)-1 before algebraic reduction:", actual_determinant_raw)
    emit("Determinant representative D(theta,x,y):", determinant_representative)
    emit("Determinant certificate H(theta,x,y):", determinant_quotient)
    emit("D(theta,x,y) - f(theta)*H(theta,x,y):", determinant_remainder)
    require(
        determinant_remainder == 0,
        "The explicit map does not have determinant identically one.",
    )
    emit("det(JF) - 1:", determinant_remainder)
    emit("det(JF) - 1 == 0:", determinant_remainder == 0)

    point_p = (sp.S.Zero, sp.S.Zero)
    point_q = (sp.S.One, sp.S.Zero)
    actual_F1_polynomial = sp.Poly(F1_actual, x, y, domain=sp.EX)
    actual_F2_polynomial = sp.Poly(F2_actual, x, y, domain=sp.EX)

    evaluations = (
        ("F1(0,0)", F1_representative, actual_F1_polynomial, point_p),
        ("F1(1,0)", F1_representative, actual_F1_polynomial, point_q),
        ("F2(0,0)", F2_representative, actual_F2_polynomial, point_p),
        ("F2(1,0)", F2_representative, actual_F2_polynomial, point_q),
    )
    formal_point_values = {}
    reduced_point_values = {}

    for label, representative, actual_polynomial, point in evaluations:
        formal_value = sp.expand(
            representative.xreplace({x: point[0], y: point[1]})
        )
        actual_raw_value = evaluate_bivariate_polynomial(
            actual_polynomial, point
        )
        require(
            sp.expand(
                actual_raw_value
                - formal_value.xreplace(root_replacement)
            ) == 0,
            f"The explicit and represented evaluations disagree for {label}.",
        )
        quotient, remainder = divide_univariate_exact(
            formal_value, theta, defining_polynomial
        )
        formal_point_values[label] = formal_value
        reduced_point_values[label] = remainder
        emit(f"{label} representative:", formal_value)
        emit(f"{label} quotient by f(theta):", quotient)
        emit(f"{label} remainder by f(theta):", remainder)
        emit(f"{label}:", sp.expand(remainder.xreplace(root_replacement)))
        require(
            remainder == 0,
            f"The required zero image fails at {label}.",
        )

    _, first_difference = divide_univariate_exact(
        formal_point_values["F1(1,0)"]
        - formal_point_values["F1(0,0)"],
        theta,
        defining_polynomial,
    )
    _, second_difference = divide_univariate_exact(
        formal_point_values["F2(1,0)"]
        - formal_point_values["F2(0,0)"],
        theta,
        defining_polynomial,
    )
    first_equal = first_difference == 0
    second_equal = second_difference == 0
    emit("F1(1,0) - F1(0,0):", first_difference)
    emit("F2(1,0) - F2(0,0):", second_difference)
    emit("F1(0,0) == F1(1,0):", first_equal)
    emit("F2(0,0) == F2(1,0):", second_equal)
    require(
        first_equal and second_equal,
        "The two distinct points do not have the same image.",
    )

    origin_substitution = {x: sp.S.Zero, y: sp.S.Zero}
    origin_jacobian_entries = []
    for derivative in (
        representative_dF1x,
        representative_dF1y,
        representative_dF2x,
        representative_dF2y,
    ):
        _, remainder = divide_univariate_exact(
            sp.expand(derivative.xreplace(origin_substitution)),
            theta,
            defining_polynomial,
        )
        origin_jacobian_entries.append(
            sp.expand(remainder.xreplace(root_replacement))
        )
    origin_jacobian = sp.Matrix(2, 2, origin_jacobian_entries)
    emit("JF(0,0):", origin_jacobian)
    emit("JF(0,0) == I_2:", origin_jacobian == sp.eye(2))
    require(
        origin_jacobian == sp.eye(2),
        "The derivative normalization failed.",
    )

    distinct_points = point_p != point_q
    emit("Collision point P:", point_p)
    emit("Collision point Q:", point_q)
    emit("P != Q:", distinct_points)
    require(distinct_points, "The purported collision points are identical.")

    emit("F(P):", (
        reduced_point_values["F1(0,0)"],
        reduced_point_values["F2(0,0)"],
    ))
    emit("F(Q):", (
        reduced_point_values["F1(1,0)"],
        reduced_point_values["F2(1,0)"],
    ))
    emit("All original coefficient equations verified exactly:", True)
    emit("All normalization, determinant, and collision checks passed:", True)
    emit("A noninjective polynomial map with det(JF) identically 1 has been certified.")


def report_unit_ideal(degree):
    emit("Groebner basis is [1].")
    emit("The ideal is the unit ideal over QQ.")
    emit(
        "No solution exists over C for the normalized Ansatz with "
        f"deg(P), deg(Q) <= {degree}."
    )
    emit("This exclusion applies only to the specified finite degree bound.")


def main(argv=None):
    degree = parse_degree(argv)
    emit("Python version:", sys.version.replace("\n", " "))
    emit("SymPy version:", sp.__version__)
    emit("Maximum nonlinear total degree d:", degree)
    emit("Arithmetic: exact rational polynomial arithmetic and exact algebraic roots.")
    emit("The two-dimensional Jacobian conjecture is open.")
    emit("Large instances may exceed available computation time or memory.")

    x, y = sp.symbols("x y")
    nonlinear_monomials = all_monomials(2, degree)
    Ppoly, Pvars = build_poly_with_vars(
        "a", x, y, nonlinear_monomials
    )
    Qpoly, Qvars = build_poly_with_vars(
        "b", x, y, nonlinear_monomials
    )
    unknowns = Pvars + Qvars
    F1 = sp.expand(x + Ppoly)
    F2 = sp.expand(y + Qpoly)

    dF1x = sp.diff(F1, x)
    dF1y = sp.diff(F1, y)
    dF2x = sp.diff(F2, x)
    dF2y = sp.diff(F2, y)
    jacobian_matrix = sp.Matrix([
        [dF1x, dF1y],
        [dF2x, dF2y],
    ])
    jacobian_determinant = sp.expand(
        dF1x * dF2y - dF1y * dF2x
    )
    determinant_minus_one = sp.expand(jacobian_determinant - 1)

    coefficient_ring = (
        sp.QQ.poly_ring(*unknowns) if unknowns else sp.QQ
    )
    determinant_polynomial = sp.Poly(
        determinant_minus_one,
        x,
        y,
        domain=coefficient_ring,
    )
    maximum_jacobian_degree = 2 * max(1, degree) - 2
    jacobian_monomials = all_monomials(0, maximum_jacobian_degree)
    jacobian_equations = []
    equation_labels = []

    for i, j in jacobian_monomials:
        coefficient = sp.expand(
            determinant_polynomial.coeff_monomial(x ** i * y ** j)
        )
        jacobian_equations.append(coefficient)
        equation_labels.append(
            f"coefficient of x^{i}*y^{j} in det(JF)-1"
        )

    point_q_substitution = {x: sp.S.One, y: sp.S.Zero}
    first_collision_equation = sp.expand(
        1 + Ppoly.xreplace(point_q_substitution)
    )
    second_collision_equation = sp.expand(
        Qpoly.xreplace(point_q_substitution)
    )
    equations = tuple(
        jacobian_equations
        + [first_collision_equation, second_collision_equation]
    )
    equation_labels.extend([
        "1 + P(1,0)",
        "Q(1,0)",
    ])

    reconstructed_determinant = sp.expand(
        sp.Add(
            *[
                equation * x ** i * y ** j
                for equation, (i, j) in zip(
                    jacobian_equations, jacobian_monomials
                )
            ]
        )
    )
    reconstruction_residual = sp.expand(
        determinant_minus_one - reconstructed_determinant
    )
    require(
        reconstruction_residual == 0,
        "The Jacobian coefficient equations omit a monomial.",
    )
    unknown_set = set(unknowns)
    for index, equation in enumerate(equations, start=1):
        require(
            equation.free_symbols.issubset(unknown_set),
            f"E_{index} contains a symbol outside the coefficient variables.",
        )

    origin_substitution = {x: sp.S.Zero, y: sp.S.Zero}
    F1_at_origin = sp.expand(F1.xreplace(origin_substitution))
    F2_at_origin = sp.expand(F2.xreplace(origin_substitution))
    jacobian_at_origin = jacobian_matrix.applyfunc(
        lambda expression: sp.expand(
            expression.xreplace(origin_substitution)
        )
    )
    require(
        F1_at_origin == 0 and F2_at_origin == 0,
        "The symbolic Ansatz does not vanish at the origin.",
    )
    require(
        jacobian_at_origin == sp.eye(2),
        "The symbolic Ansatz does not have identity derivative at the origin.",
    )

    emit("Spatial variables:", (x, y))
    emit(
        "Nonlinear monomials:",
        [x ** i * y ** j for i, j in nonlinear_monomials],
    )
    for coefficient, (i, j) in zip(Pvars, nonlinear_monomials):
        emit("P coefficient:", coefficient, "monomial:", x ** i * y ** j)
    for coefficient, (i, j) in zip(Qvars, nonlinear_monomials):
        emit("Q coefficient:", coefficient, "monomial:", x ** i * y ** j)

    emit("Unknown coefficients in lexicographic priority order:", unknowns)
    emit("Number of unknown coefficients:", len(unknowns))
    emit("P(x,y):", Ppoly)
    emit("Q(x,y):", Qpoly)
    emit("F1(x,y):", F1)
    emit("F2(x,y):", F2)
    emit("JF(x,y):", jacobian_matrix)
    emit("Expanded det(JF):", jacobian_determinant)
    emit("Expanded det(JF)-1:", determinant_minus_one)
    emit("Jacobian coefficient reconstruction residual:", reconstruction_residual)
    emit("F1(0,0):", F1_at_origin)
    emit("F2(0,0):", F2_at_origin)
    emit("JF(0,0):", jacobian_at_origin)
    emit("Number of equations, including identically zero equations:", len(equations))
    emit("Number of nonzero equations:", sum(equation != 0 for equation in equations))

    for index, (label, equation) in enumerate(
        zip(equation_labels, equations),
        start=1,
    ):
        emit(f"E_{index} [{label}]:", equation, "= 0")

    if not unknowns:
        certify_constant_system(equations)
        report_unit_ideal(degree)
        return 0

    emit("Computing the reduced lexicographic Groebner basis over QQ.")
    basis = rational_groebner(equations, unknowns)
    print_basis("Reduced Groebner basis G:", basis)
    certify_reduced_groebner(basis, equations)

    if is_unit_basis(basis):
        report_unit_ideal(degree)
        return 0

    emit("The certified ideal is proper over QQ.")
    emit("Constructing one exact algebraic coefficient assignment.")
    theta, defining_polynomial, alpha, representatives = (
        construct_exact_algebraic_point(basis)
    )
    verify_and_print_solution(
        degree,
        x,
        y,
        F1,
        F2,
        unknowns,
        equations,
        theta,
        defining_polynomial,
        alpha,
        representatives,
    )
    return 0


def entrypoint():
    try:
        if hasattr(sys, "set_int_max_str_digits"):
            sys.set_int_max_str_digits(0)
        return main()
    except KeyboardInterrupt:
        print(
            "\nINCOMPLETE: computation interrupted. "
            "Interruption is not evidence of nonexistence.",
            file=sys.stderr,
            flush=True,
        )
        return 130
    except MemoryError:
        print(
            "INCOMPLETE: available memory was exhausted before the workflow finished. "
            "Resource exhaustion is not evidence of nonexistence.",
            file=sys.stderr,
            flush=True,
        )
        return 3
    except BrokenPipeError:
        sys.stdout = open(os.devnull, "w", encoding="utf-8")
        return 1
    except Exception:
        print(
            "FAILED: the requested computation or exact certification did not "
            "complete normally. No unverified candidate is certified.",
            file=sys.stderr,
            flush=True,
        )
        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(entrypoint())
