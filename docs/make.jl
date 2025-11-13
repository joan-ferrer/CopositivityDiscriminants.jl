using CopositivityDiscriminants
using Documenter

DocMeta.setdocmeta!(CopositivityDiscriminants, :DocTestSetup, :(using CopositivityDiscriminants); recursive=true)

makedocs(;
    modules=[CopositivityDiscriminants],
    authors="Elisenda Feliu, Joan Ferrer and Máté Telek",
    sitename="CopositivityDiscriminants.jl",
    format=Documenter.HTML(;
        canonical="https://joan-ferrer.github.io/CopositivityDiscriminants.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/joan-ferrer/CopositivityDiscriminants.jl",
    devbranch="main",
)
