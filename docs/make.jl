using Documenter

push!(LOAD_PATH, joinpath(@__DIR__,  "../src/") )
using MultiObjectiveProblems

makedocs(sitename="MultiObjectiveProblems Documentation")