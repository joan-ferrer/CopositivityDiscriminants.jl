using HomotopyContinuation
using Oscar
using Arblib
using .CopositivityDiscriminants: CoposCheckResult  

# helper to get certified interval of a certificate
get_interval(c) = begin
    if c isa HomotopyContinuation.ExtendedSolutionCertificate
        X = HomotopyContinuation.certified_solution_interval_after_krawczyk(c)
        X === nothing && (X = HomotopyContinuation.certified_solution_interval(c))
        return X
    else
        return HomotopyContinuation.certified_solution_interval(c)
    end
end

"""
    check_copositivity(f::Expression, nonseparable::Bool=false)

Builds f_t by multiplying all negative-coefficient monomials by a parameter t,
sets up the system [f_t; x_i * ∂f_t/∂x_i], solves & certifies, and
returns a `CoposCheckResult` summarizing:
- whether the minimum positive solution has t ≥ 1,
- whether t=1 lies in that min-t certificate's certified interval,
- the min-t estimate, certificates, and system.

`nonseparable=true` constructs an special homotopy for polynomials with 
nonseparable signed supports: it constructs barycentric start parameters via an Oscar polyhedron,
tracks to the target coefficients, substitutes parameters, and certifies.
"""
function check_copositivity(f::Expression, nonseparable::Bool=false)
    # Variables / data from f
    vars = variables(f)
    dimension = length(vars)

    # HC variable t first, then the original vars
    HomotopyContinuation.@var t
    vars_t = [t; vars...]

    # Exponents & coefficients
    exps, coeffs = exponents_coefficients(f, vars)   # exps: (dimension × nterms)


    # ---------------- Oscar checks on the Newton polytope ----------------

    # Newton polytope: convex hull of all exponent vectors
    # exps has columns = exponent vectors, so we transpose to get rows = points
    newton_poly = Oscar.convex_hull(Matrix(exps'))

    # 1) Hull must be full-dimensional
    @assert Oscar.is_fulldimensional(newton_poly) "Newton polytope of f must be full-dimensional."

    # Split by sign 
    pos_idx = findall(>=(0), coeffs)
    neg_idx = findall(<(0),  coeffs)

    # 2) Every negative exponent must be an interior lattice point

    # Check that there is at least one negative term
    @assert !isempty(neg_idx) "Polynomial must have at least one term with negative coefficient."

    interior_pts = collect(Oscar.interior_lattice_points(newton_poly))

    @assert !isempty(interior_pts) "Newton polytope has no interior lattice points but f has monomials with negative coefficients."

    for j in neg_idx
        e = exps[:, j]                         # exponent vector (Int)
        e_zz = [Oscar.ZZ(e[i]) for i in 1:dimension]  # same vector in ZZ

        is_interior = any(p -> all(p[i] == e_zz[i] for i in 1:dimension),
                        interior_pts)

        @assert is_interior "Exponent $(collect(e)) with negative coefficient is not an interior lattice point of the Newton polytope."
    end


    pos_coeffs = coeffs[pos_idx]
    neg_coeffs = coeffs[neg_idx]

    # Collect exponent columns as small vectors for easy reuse
    pos_exps = [exps[:, j] for j in pos_idx]
    neg_exps = [exps[:, j] for j in neg_idx]

    # Helpers to build monomials quickly
    monomial(vs, e) = prod(vs[i]^e[i] for i in 1:length(vs))
    pos_monos = [monomial(vars, e) for e in pos_exps]
    neg_monos = [monomial(vars, e) for e in neg_exps]

    if !nonseparable
        # ------------------------ General case ------------------------
        # f_t = sum(pos) + t * sum(neg)
        f_t = sum(pos_coeffs[j] * pos_monos[j] for j in eachindex(pos_monos)) +
              t * sum(neg_coeffs[j] * neg_monos[j] for j in eachindex(neg_monos))

    
        eqs = [f_t; [v * differentiate(f_t, v) for v in vars]...]

        F = System(eqs; variables=vars_t)
        result = HomotopyContinuation.solve(F)

        C = certify(F, result)
        certs = certificates(C)

        # Keep only positive certificates
        pos = [c for c in certs if HomotopyContinuation.is_positive(c)]

        # If no positive certificates
        if isempty(pos)
            return CoposCheckResult(;
                copositive = false,
                method = :general,
                t_min = NaN,
                system = F,
                cert_result = C,
                positive_certs = HomotopyContinuation.AbstractSolutionCertificate[],
                certified_interval_t_min = nothing,
            )
        end

        # Choose the positive cert with minimal real midpoint of t (first coord)
        t_mids = map(pos) do c
            z = solution_candidate(c)      # Complex midpoint vector
            real(z[1])
        end
        imin = argmin(t_mids)
        cmin = pos[imin]
        t_min = t_mids[imin]

        X = get_interval(cmin)
        tball = X === nothing ? nothing : X[1]

        has_t1 = false
        if tball !== nothing
            has_t1 = Arblib.contains(Arblib.real(tball), 1) &&
                     Arblib.contains(Arblib.imag(tball), 0)
        end

        if has_t1
            return CoposCheckResult(;
            copositive = missing,
            method = :general,
            t_min = t_min,
            system = F,
            cert_result = C,
            positive_certs = pos,
            certified_interval_t_min = tball,
        )
        else
            return CoposCheckResult(;
                copositive = t_min ≥ 1,
                method = :general,
                t_min = t_min,
                system = F,
                cert_result = C,
                positive_certs = pos,
                certified_interval_t_min = tball,
            )
        end        

    else
        # ----------------  Nonseparable case ----------------

        num_pos = length(pos_exps)
        num_neg = length(neg_exps)

        @assert num_neg ≥ 1 "Nonseparable=true assumes at least one negative term."

        Epos = (num_pos == 0) ? zeros(Int, dimension, 0) : hcat(pos_exps...)
        col_pos = vcat(ones(Int, 1, num_pos), Epos)          # (dimension+1)×num_pos
        col_neg = -vcat(1, neg_exps[1])                      # (dimension+1)
        A = hcat(col_pos, col_neg)                           # (dimension+1)×(num_pos+1)

        # Polyhedron: P = {y ≥ 0, A*y = 0}
        Iblock = Matrix{Int}(I, num_pos + 1, num_pos + 1)
        P = Oscar.polyhedron((-Iblock, zeros(Int, num_pos + 1)), (A, zeros(Int, dimension + 1)))
        interior = Oscar.relative_interior_point(P)          # rational vector of length num_pos+1

        # Convert rationals -> Float64 
        vec_numer = BigInt[numerator(v) for v in interior]
        vec_denom = BigInt[denominator(v) for v in interior]
        interior_float = vec_numer ./ vec_denom              # barycentric coordinates

        # Initial parameter coefficients p for the start system:
        #   - positive terms take barycentric weights,
        #   - ONLY the first negative term is active with weight -last_bary.
        total_terms = length(coeffs)
        initial_coeffs = fill(0.0, total_terms)
        for j in 1:num_pos
            initial_coeffs[j] = interior_float[j]
        end
        initial_coeffs[num_pos + 1] = -interior_float[end]

        # Parameterized homotopy in the coefficients p[1:total_terms]
        HomotopyContinuation.@var p[1:total_terms]

        f_t = sum(p[j] * pos_monos[j] for j in 1:num_pos) +
              t * sum(p[num_pos + j] * neg_monos[j] for j in 1:num_neg)

        eqs = [f_t; [v * differentiate(f_t, v) for v in vars]...]
        F = System(eqs; variables=vars_t, parameters=p)

        # Start & target
        start_solutions = ones(Float64, dimension + 1)
        reordered_coeffs = vcat(pos_coeffs, neg_coeffs)

        res = HomotopyContinuation.solve(
            F, [start_solutions];
            start_parameters = Float64.(initial_coeffs),
            target_parameters = reordered_coeffs,
            seed = 0x68a5c2c6,
        )

        # Substitute target parameters and certify
        eqs_sub = subs(eqs, Dict(p .=> reordered_coeffs))
        F_sub = HomotopyContinuation.System(eqs_sub)

        C = certify(F_sub, solutions(res))
        certs = certificates(C)
        pos = [c for c in certs if HomotopyContinuation.is_positive(c)]

        if isempty(pos)
            return CoposCheckResult(;
                copositive = false,
                method = :nonseparable,
                t_min = NaN,
                system = F_sub,
                cert_result = C,
                positive_certs = HomotopyContinuation.AbstractSolutionCertificate[],
                certified_interval_t_min= nothing,
            )
        end

        # Min by real midpoint of t among positive certificates
        t_mids = map(pos) do c
            z = solution_candidate(c)
            real(z[1])
        end
        imin = argmin(t_mids)
        cmin = pos[imin]
        t_min = t_mids[imin]

        X = get_interval(cmin)
        tball = X === nothing ? nothing : X[1]

        has_t1 = false
        if tball !== nothing
            has_t1 = Arblib.contains(Arblib.real(tball), 1) &&
                     Arblib.contains(Arblib.imag(tball), 0)
        end

        if has_t1
            return CoposCheckResult(;
            copositive = missing,
            method = :nonseparable,
            t_min = t_min,
            system = F,
            cert_result = C,
            positive_certs = pos,
            certified_interval_t_min = tball,
        )
        else
            return CoposCheckResult(;
                copositive = t_min ≥ 1,
                method = :nonseparable,
                t_min = t_min,
                system = F,
                cert_result = C,
                positive_certs = pos,
                certified_interval_t_min = tball,
            )
        end
    end
end
