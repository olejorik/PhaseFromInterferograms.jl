using PhaseFromInterferograms
using Documenter, Literate

DocMeta.setdocmeta!(
    PhaseFromInterferograms, :DocTestSetup, :(using PhaseFromInterferograms); recursive=true
)

# TODO Generate tutorials from literate files
@info "current dir =$(@__DIR__)"
tutorials_folder = (@__DIR__) * "/../tutorials"
docs_tutorials_folder = (@__DIR__) * "/src/tutorials"
@info tutorials_folder
for f in readdir(tutorials_folder; join=true)
    Literate.markdown(f, docs_tutorials_folder)
end


makedocs(;
    modules=[PhaseFromInterferograms],
    authors="Oleg Soloviev",
    repo="https://github.com/olejorik/PhaseFromInterferograms.jl/blob/{commit}{path}#L{line}",
    sitename="PhaseFromInterferograms.jl",
    # doctest=:fix,
    checkdocs=:exports,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://olejorik.github.io/PhaseFromInterferograms.jl",
        assets=String[],
    ),
    clean=false,
    pages=[
        "Home" => "index.md", "Manual" => ["Finding tilts" => "tutorials/FindingTilts.md"]
    ],
)

deploydocs(; repo="github.com/olejorik/PhaseFromInterferograms.jl")
