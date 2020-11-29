using MultiObjectiveProblems
using Documenter

DocMeta.setdocmeta!(MultiObjectiveProblems, :DocTestSetup, :(using MultiObjectiveProblems); recursive=true)

makedocs(;
    modules=[MultiObjectiveProblems],
    authors="Manuel Berkemeier <manuelbb@math.upb.de>",
    repo="https://git.cs.uni-paderborn.de/manuelbb/MultiObjectiveProblems.jl/blob/{commit}{path}#{line}",
    sitename="MultiObjectiveProblems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://manuelbb.gitlab.io/MultiObjectiveProblems.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
