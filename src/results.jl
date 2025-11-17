using Base: @kwdef

const MaybeBool = Union{Bool,Missing}

"""
    CoposCheckResult

Structured summary of the result returned by [`check_copositivity`](@ref).

Instances are constructed via keyword arguments and
contain both the numerical estimate of the minimal parameter `t` and the
certificates produced by HomotopyContinuation.

# Fields
- `copositive::MaybeBool`  
  Tri-state outcome of the copositivity check:
  - `true`  – the procedure certifies that the input is copositive;
  - `false` – the procedure certifies that it is **not** copositive;
  - `missing` – the procedure was inconclusive (e.g. certification failed).

- `method::Symbol`  
  Method used for the check:
  - `:general`      – Computes all solution of the system and minimal t-value among the positive ones.
  - `:nonseparable` – specialised homotopy for nonseparable signed supports.

- `t_min::Float64`  
  Midpoint estimate of the minimal parameter `t` found by the homotopy
  continuation setup (the “min-t” value referred to in the documentation of
  [`check_copositivity`](@ref)).

- `system::Union{Nothing, HomotopyContinuation.System}`  
  The `HomotopyContinuation.System` that was solved in the copositivity check,
  or `nothing` if it was not stored.

- `cert_result::Union{Nothing, HomotopyContinuation.CertificationResult}`  
  Overall certification result returned by HomotopyContinuation, or `nothing`
  if certification was not attempted or not kept.

- `positive_certs::Vector{HomotopyContinuation.AbstractSolutionCertificate}`  
  Certificates for the positive solutions that were used to decide
  copositivity. This vector may be empty if no certified solutions were
  obtained.

- `certified_interval_t_min::Any`  
  Rigorous enclosure of the true minimal `t` (typically an interval-like
  object), or `nothing` if no certified interval was produced.

The convenience alias `MaybeBool = Union{Bool,Missing}` is used for the
`copositive` field.
"""
@kwdef struct CoposCheckResult
    copositive::MaybeBool                      # true | false | missing
    method::Symbol                         # :general | :nonseparable
    t_min::Float64                         # midpoint estimate of min t
    system::Union{Nothing,HomotopyContinuation.System} = nothing # system used for check
    cert_result::Union{Nothing,HomotopyContinuation.CertificationResult} = nothing
    positive_certs::Vector{HomotopyContinuation.AbstractSolutionCertificate} =
        HomotopyContinuation.AbstractSolutionCertificate[]
    certified_interval_t_min::Any = nothing          
end