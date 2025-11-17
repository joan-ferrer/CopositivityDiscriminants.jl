```@meta
CurrentModule = CopositivityDiscriminants
```

# CopositivityDiscriminants

[CopositivityDiscriminants](https://github.com/joan-ferrer/CopositivityDiscriminants.jl) is a proof-of-concept Julia package that implements the methods developed in LINK TO PAPER for checking copositivity (i.e. nonnegativity on  the positive orthant) of real polynomials. 

Currently, it implements those methods for polynomials with full dimensional Newton polytope and such that no term with negative coefficient lies on the boundary of the Newton polytope. Using the notation from LINK TO PAPER, it supports polynomials $f\in\mathbb{R}[x_1,\dots,x_n]$ whose signed support $(\mathcal{A}_+,\mathcal{A}_-)$ satisfies $dim(conv(\mathcal{A}_+))=n$ and $\mathcal{A}_-\subseteq int(conv(\mathcal{A}_+))$.



```@index
```

## Installation
Use the following Julia commands to intall the package.

```julia
using Pkg
Pkg.add(url="https://github.com/joan-ferrer/CopositivityDiscriminants.jl.git")
```


```@autodocs
Modules = [CopositivityDiscriminants]
```



