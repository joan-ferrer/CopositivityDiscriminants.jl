# CopositivityDiscriminants

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://joan-ferrer.github.io/CopositivityDiscriminants.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://joan-ferrer.github.io/CopositivityDiscriminants.jl/dev/)
[![Build Status](https://github.com/joan-ferrer/CopositivityDiscriminants.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/joan-ferrer/CopositivityDiscriminants.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/joan-ferrer/CopositivityDiscriminants.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/joan-ferrer/CopositivityDiscriminants.jl)


[CopositivityDiscriminants.jl](https://github.com/joan-ferrer/CopositivityDiscriminants.jl) is a proof-of-concept Julia package that implements the methods developed in LINK TO PAPER for checking copositivity (i.e. nonnegativity on  the positive orthant) of real polynomials. 

Currently, it implements those methods for polynomials with full dimensional Newton polytope and such that no term with negative coefficient lies on the boundary of the Newton polytope. Using the notation from LINK TO PAPER, it supports polynomials $f\in\mathbb{R}[x_1,\dots,x_n]$ whose signed support $ (\mathcal{A}_+,\mathcal{A}_-) $ satisfies $ dim(conv(\mathcal{A}_+))=n$ and $\mathcal{A}_-\subseteq int(conv(\mathcal{A}_+))$.

## Installation
Use the following Julia commands to intall the package.

```julia
using Pkg
Pkg.add(url="https://github.com/joan-ferrer/CopositivityDiscriminants.jl.git")
```