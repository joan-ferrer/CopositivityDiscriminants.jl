using LinearAlgebra
using Oscar
import HomotopyContinuation as HC

# --- support from HC (strict sign split) ---
function _split_support_from_hc_expression(p)
    vars = HC.variables(p)
    M, coeffs = HC.ModelKit.exponents_coefficients(p, vars; expanded=false)  
    E = permutedims(Matrix(M))                                                
    pos_ix = findall(>(0), coeffs)
    neg_ix = findall(<(0), coeffs)
    return E[pos_ix, :], E[neg_ix, :]
end

# --- small helpers ---
@inline _barycentric(Ainv::AbstractMatrix, x::AbstractVector) = Ainv * vcat(Float64.(x), 1.0)
@inline _in_simplex_nonstrict(λ; tol=1e-9) = all(λ .>= -tol)

# Convert any “point-like” vector to a tuple of Int for hashable set membership
_to_int_tuple(v) = Tuple(Int.(collect(v)))

function _negatives_subset_of_interior_oscar(Apos::AbstractMatrix, Aneg::AbstractMatrix)::Bool
    mpos = size(Apos, 1)

    # Build the convex hull 
    pts = [Oscar.point_vector(Oscar.QQ, vec(Apos[i, :])) for i in 1:mpos]
    P   = Oscar.convex_hull(pts)

    # Collect all interior lattice points as tuples of Int (any length)
    ILPset = Set{Tuple{Vararg{Int}}}()
    for q in Oscar.interior_lattice_points(P)
        # 'q' is iterable over coordinates (often Nemo fmpz); convert to Ints
        push!(ILPset, Tuple(Int.(collect(q))))
    end

    # Every negative exponent must be one of those interior lattice points
    for i in 1:size(Aneg, 1)
        x = Tuple(Int.(vec(Aneg[i, :])))  # normalize Int32/Int64 → Int
        if !(x in ILPset)
            return false
        end
    end
    return true
end

# Core check across all regular triangulations 
function _nonsep_core_all_regular_oscar(Apos::AbstractMatrix, Aneg::AbstractMatrix;
        tol::Real=1e-9, verbose::Bool=true)::Bool

    isempty(Apos) && return false
    n = size(Apos, 2)
    size(Apos,1) < n+1 && return false

    # 1) full-dimensional hull?
    pts_hull = [Oscar.point_vector(Oscar.QQ, vec(Apos[i,:])) for i in 1:size(Apos,1)]
    P = Oscar.convex_hull(pts_hull)
    if Oscar.dim(P) < Oscar.ambient_dim(P)
        verbose && @warn "Convex hull of positive exponents is not full dimensional."
        return false
    end

    # 2) negatives strictly inside via interior_lattice_points
    _negatives_subset_of_interior_oscar(Apos, Aneg) || begin
        verbose && @info "The signed support has negative terms on the boundary"
        return false
    end

    # 3) per-simplex condition for ALL regular triangulations
    pts  = [Oscar.point_vector(Oscar.QQ, vec(Apos[i,:])) for i in 1:size(Apos,1)]
    tris = Oscar.regular_triangulations(pts; full=false)   # 1-based simplices
    k    = n + 1
    Tneg = size(Aneg,1)

    for cells in tris                 # each triangulation
        for idxs in cells             # each simplex (indices into Apos)
            length(idxs) == k || continue
            V = transpose(Apos[idxs, :]) |> x -> Float64.(x)
            A = vcat(V, ones(1, k))
            rank(A) == k || continue
            Ainv = inv(A)

            # negatives contained (allow boundary)
            neg_in = Int[]
            lambdas = Vector{Vector{Float64}}()
            for i in 1:Tneg
                λ = _barycentric(Ainv, vec(Aneg[i,:]))
                if _in_simplex_nonstrict(λ; tol=tol)
                    push!(neg_in, i); push!(lambdas, λ)
                end
            end

            isempty(neg_in) && continue           # case (1)
            if length(neg_in) == Tneg; continue; end  # case (2a)

            # case (2b): all negatives here lie on the SAME facet
            zero_sets = [Set(findall(j -> abs(λ[j]) ≤ tol, 1:k)) for λ in lambdas]
            isempty(reduce(intersect, zero_sets)) && return false
        end
    end
    return true
end

# --- public API ---

"""
    nonseparable_support(p; tol::Real = 1e-9, verbose::Bool = true) -> Bool

Given a polynomial `p` (a `HomotopyContinuation.Expression`), compute the signed
support of `p` and decide whether it is *nonseparable* in the sense used in
this package.

The polynomial is written as
\\[
p(x) = \\sum_{\\alpha \\in \\mathcal{A}_{+}} c_{\\alpha} x^{\\alpha}
      + \\sum_{\\beta \\in \\mathcal{A}_{-}} c_{\\beta} x^{\\beta},
\\]
where `Apos`(\\(\\mathcal{A}_+\\)) are the exponents with positive coefficients and `Aneg` (\\(\\mathcal{A}_-\\) ) the
exponents with negative coefficients. Internally

- `_split_support_from_hc_expression(p)` builds the exponent matrices
  `Apos` and `Aneg`,
- Oscar is used to form the convex hull of the rows of `Apos`,
- negative exponents are required to be *strictly interior* lattice points of
  this hull, and
- for **every** regular triangulation of the positive support, each simplex
  either
  1. contains no negative exponent,
  2. contains all negative exponents, or
  3. contains some but not all negative exponents, and in this case all
     such negatives lie on a common facet (detected via barycentric
     coordinates).

If all these conditions are satisfied, the signed support is declared
nonseparable.

# Arguments
- `p`: a `HomotopyContinuation.Expression` whose signed support is to be
  tested.
- `tol`: numerical tolerance used when testing barycentric coordinates for
  simplex membership and detecting whether a point lies on a facet.
- `verbose`: if `true`, prints warnings / info messages when the positive
  support is not full dimensional or when negative exponents lie on the
  boundary.

# Returns
- `true` if the signed support of `p` is nonseparable according to the above
  criteria;
- `false` otherwise (including degenerate cases such as non–full-dimensional
  positive support or negative exponents not strictly interior).
"""

function nonseparable_support(p; tol::Real=1e-9, verbose::Bool=true)::Bool
     Apos, Aneg = _split_support_from_hc_expression(p)
    _nonsep_core_all_regular_oscar(Apos, Aneg; tol=tol, verbose=verbose)
end
