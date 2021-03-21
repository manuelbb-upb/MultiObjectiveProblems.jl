using Pkg
Pkg.activate(@__DIR__)
#%%
using MultiObjectiveProblems
using Documenter

DocMeta.setdocmeta!(MultiObjectiveProblems, :DocTestSetup, :(using MultiObjectiveProblems); recursive=true)

makedocs(;
    modules=[MultiObjectiveProblems],
    authors="Manuel Berkemeier",
    repo="https://github.com/manuelbb-upb/MultiObjectiveProblems.jl/blob/{commit}{path}#{line}",
    sitename="MultiObjectiveProblems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://manuelbb-upb.github.io/MultiObjectiveProblems.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "DTLZ Problems" => "dtlz.md",
        "ZDT Problems" => "zdt.md",
        "MHT Problems" => "mht.md",
    ],
)

deploydocs(;
    repo="github.com/manuelbb-upb/MultiObjectiveProblems.jl",
)
