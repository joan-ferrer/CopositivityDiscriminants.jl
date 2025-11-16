using LinearAlgebra
using Oscar
import HomotopyContinuation as HC

# --- support from HC (strict sign split; no tolerance) ---
function _split_support_from_hc_expression(p)
    vars = HC.variables(p)
    M, coeffs = HC.ModelKit.exponents_coefficients(p, vars; expanded=false)  # n×m
    E = permutedims(Matrix(M))                                                # m×n
    pos_ix = findall(>(0), coeffs)
    neg_ix = findall(<(0), coeffs)
    return E[pos_ix, :], E[neg_ix, :]
end

# --- small helpers ---
@inline _barycentric(Ainv::AbstractMatrix, x::AbstractVector) = Ainv * vcat(Float64.(x), 1.0)
@inline _in_simplex_nonstrict(λ; tol=1e-9) = all(λ .>= -tol)

# Convert any “point-like” vector to a tuple of Int for hashable set membership
_to_int_tuple(v) = Tuple(Int.(collect(v)))

function _negatives_subset_of_interior_oscar(Ppos::AbstractMatrix, Pneg::AbstractMatrix)::Bool
    mpos = size(Ppos, 1)

    # Build the hull once
    pts = [Oscar.point_vector(Oscar.QQ, vec(Ppos[i, :])) for i in 1:mpos]
    P   = Oscar.convex_hull(pts)

    # Collect all interior lattice points as tuples of Int (any length)
    ILPset = Set{Tuple{Vararg{Int}}}()
    for q in Oscar.interior_lattice_points(P)
        # 'q' is iterable over coordinates (often Nemo fmpz); convert robustly to Ints
        push!(ILPset, Tuple(Int.(collect(q))))
    end

    # Every negative exponent must be one of those interior lattice points
    for i in 1:size(Pneg, 1)
        x = Tuple(Int.(vec(Pneg[i, :])))  # normalize Int32/Int64 → Int
        if !(x in ILPset)
            return false
        end
    end
    return true
end

# Core check across ALL regular triangulations (OSCAR only)
function _nonsep_core_all_regular_oscar(Ppos::AbstractMatrix, Pneg::AbstractMatrix;
        tol::Real=1e-9, verbose::Bool=true)::Bool

    isempty(Ppos) && return false
    n = size(Ppos, 2)
    size(Ppos,1) < n+1 && return false

    # 1) full-dimensional hull?
    pts_hull = [Oscar.point_vector(Oscar.QQ, vec(Ppos[i,:])) for i in 1:size(Ppos,1)]
    P = Oscar.convex_hull(pts_hull)
    if Oscar.dim(P) < Oscar.ambient_dim(P)
        verbose && @warn "Convex hull of positive exponents is not full dimensional."
        return false
    end

    # 2) negatives strictly inside via interior_lattice_points
    _negatives_subset_of_interior_oscar(Ppos, Pneg) || begin
        verbose && @info "The signed support is not full dimensional. You may want to use the function convert_support_to_full_dimension()."
        return false
    end

    # 3) per-simplex condition for ALL regular triangulations (same as before)
    pts  = [Oscar.point_vector(Oscar.QQ, vec(Ppos[i,:])) for i in 1:size(Ppos,1)]
    tris = Oscar.regular_triangulations(pts; full=false)   # 1-based simplices
    k    = n + 1
    Tneg = size(Pneg,1)

    for cells in tris                 # each triangulation
        for idxs in cells             # each simplex (indices into Ppos)
            length(idxs) == k || continue
            V = transpose(Ppos[idxs, :]) |> x -> Float64.(x)
            A = vcat(V, ones(1, k))
            rank(A) == k || continue
            Ainv = inv(A)

            # negatives contained (allow boundary)
            neg_in = Int[]
            lambdas = Vector{Vector{Float64}}()
            for i in 1:Tneg
                λ = _barycentric(Ainv, vec(Pneg[i,:]))
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
function nonseparable_support(p; tol::Real=1e-9, verbose::Bool=true)::Bool
    Ppos, Pneg = _split_support_from_hc_expression(p)
    _nonsep_core_all_regular_oscar(Ppos, Pneg; tol=tol, verbose=verbose)
end
