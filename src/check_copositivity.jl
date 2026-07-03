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
    check_copositivity(f::Expression, nonseparable::Bool=false,h::Union{Nothing, AbstractVector{Int}}=nothing,use_extended_cert::Bool=false)

Builds f_t by multiplying all negative-coefficient monomials by a parameter t,
sets up the system [f_t; x_i * ∂f_t/∂x_i], solves & certifies, and
returns a `CoposCheckResult` summarizing:
- whether the minimum positive solution has t ≥ 1,
- whether t=1 lies in that min-t certificate's certified interval,
- the min-t estimate, certificates, and system.

Keyword Arguments:
-`nonseparable::Bool` constructs an special homotopy for polynomials with 
nonseparable signed supports: it constructs barycentric start parameters via an Oscar polyhedron,
tracks to the target coefficients, substitutes parameters, and certifies. It is false by default.
-`h::Union{Nothing, AbstractVector{Int}}`: Specifies the the height function determining the exponent
 of t multiplying each negative term. The order should match the order of the negative terms determined 
 by HomotopyContinuation's internal polynomial parsing (typically graded lexicographic). 
 If `nothing`, defaults to 1 for all negative terms.
 - `use_extended_cert::Bool`: If true, uses extended certificates to try to further shrink the certification interval. Default is false.
"""
function check_copositivity(f::Expression; nonseparable::Bool=false, h::Union{Nothing, AbstractVector{Int}}=nothing)
    vars = variables(f)
    dimension = length(vars)

    # HC variable t first, then the original vars
    HomotopyContinuation.@var t
    vars_t = [t; vars...]

    exps, coeffs = exponents_coefficients(f, vars)   # exps: (dimension × nterms)


    # ---------------- Oscar checks on the Newton polytope ----------------

    # exps has columns = exponent vectors, so we transpose to get rows = points
    newton_poly = Oscar.convex_hull(Matrix(exps'))

    # 1) Hull must be full-dimensional
    @assert Oscar.is_fulldimensional(newton_poly) "Newton polytope of f must be full-dimensional."

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


    # ---------------- Handle t exponents ----------------
    if h === nothing
        t_exps = ones(Int, length(neg_idx))
    else
        @assert length(h) == length(neg_idx) "Provided $(length(h)) h exponents, but f has $(length(neg_idx)) negative terms."
        t_exps = h

        mapping_strings = ["  t^$(t_exps[j]) * ($(neg_coeffs[j]) * $(neg_monos[j]))" for j in eachindex(neg_monos)]
        info_message = "Applied height function to negative terms:\n" * join(mapping_strings, "\n")
        @info info_message
    end
    local final_system, C, pos_certs, method_used
    if !nonseparable
        # ------------------------ General case ------------------------
        method_used = :general
        
        # Apply the specific t exponent to each negative term
        f_t = sum(pos_coeffs[j] * pos_monos[j] for j in eachindex(pos_monos)) +
              sum((t^t_exps[j]) * neg_coeffs[j] * neg_monos[j] for j in eachindex(neg_monos))

        eqs = [f_t; [v * differentiate(f_t, v) for v in vars]...]
        final_system = System(eqs; variables=vars_t)
        
        result = HomotopyContinuation.solve(final_system)
        C = certify(final_system, result,extended_certificates=use_extended_cert)
        
    else
        # ----------------  Nonseparable case ----------------
        method_used = :nonseparable
        num_pos = length(pos_exps)
        num_neg = length(neg_exps)

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
              sum((t^t_exps[j]) * p[num_pos + j] * neg_monos[j] for j in 1:num_neg)

        eqs = [f_t; [v * differentiate(f_t, v) for v in vars]...]
        F = System(eqs; variables=vars_t, parameters=p)

        start_solutions = ones(Float64, dimension + 1)
        reordered_coeffs = vcat(pos_coeffs, neg_coeffs)

        res = HomotopyContinuation.solve(
            F, [start_solutions];
            start_parameters = Float64.(initial_coeffs),
            target_parameters = reordered_coeffs,
            seed = 0x68a5c2c6,
        )

        eqs_sub = subs(eqs, Dict(p .=> reordered_coeffs))
        final_system = HomotopyContinuation.System(eqs_sub) 
        C = certify(final_system, solutions(res),extended_certificates=use_extended_cert)
    end
# ----------------  Extraction and Return ----------------
    certs = certificates(C)
    pos_certs = [c for c in certs if HomotopyContinuation.is_positive(c)]

    if isempty(pos_certs)
        return CoposCheckResult(;
            copositive = false,
            method = method_used,
            system = final_system,
            cert_result = C,
            positive_certs = HomotopyContinuation.AbstractSolutionCertificate[],
            t_min_interval = nothing,
        )
    end

    intervals = map(get_interval, pos_certs)
    t_mids = map(intervals) do X
        X === nothing ? Inf : Arblib.midref(Arblib.real(X[1]))
    end
    
    imin = argmin(t_mids)
    cmin = pos_certs[imin]
    
    X_min = intervals[imin]
    tball = X_min === nothing ? nothing : X_min[1]

    has_t1 = false
    is_strictly_greater = false

    if tball !== nothing
        real_part = Arblib.real(tball)
        
        has_t1 = Arblib.contains(real_part, 1) && Arblib.contains(Arblib.imag(tball), 0)
        
        if !has_t1
            is_strictly_greater = real_part > 1
        end
    end

    return CoposCheckResult(;
        copositive = has_t1 ? missing : is_strictly_greater,
        method = method_used,
        system = final_system,
        cert_result = C,
        positive_certs = pos_certs,
        t_min_interval = tball,
    )
end