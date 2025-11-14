using Base: @kwdef

const MaybeBool = Union{Bool,Missing}

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