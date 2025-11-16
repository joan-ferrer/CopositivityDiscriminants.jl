module CopositivityDiscriminants

using HomotopyContinuation
using Oscar
using LinearAlgebra
using Arblib

export check_copositivity, CoposCheckResult, nonseparable_support

include("results.jl")
include("check_copositivity.jl")
include("nonseparable_support.jl")

end
