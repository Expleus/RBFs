using RBFs
using Documenter

DocMeta.setdocmeta!(RBFs, :DocTestSetup, :(using RBFs); recursive=true)

makedocs(;
    modules=[RBFs],
    authors="EIA-FTA",
    sitename="RBFs.jl",
    format=Documenter.HTML(;
        canonical="https://Expleus.github.io/RBFs.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Expleus/RBFs.jl",
    devbranch="master",
)
