using ALGCM1D
using Documenter

DocMeta.setdocmeta!(ALGCM1D, :DocTestSetup, :(using ALGCM1D); recursive=true)

makedocs(;
    modules=[ALGCM1D],
    authors="Udi Strobach <udistr@gmail.com> and Ziv Moreno <moreno_ziv@yahoo.com>",
    repo="https://github.com/udistr/ALGCM1D.jl/blob/{commit}{path}#{line}",
    sitename="ALGCM1D.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
