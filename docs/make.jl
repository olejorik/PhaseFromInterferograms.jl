using PhaseFromInterferograms
using Documenter

makedocs(;
    modules=[PhaseFromInterferograms],
    authors="Oleg Soloviev",
    repo="https://github.com/olejorik/PhaseFromInterferograms.jl/blob/{commit}{path}#L{line}",
    sitename="PhaseFromInterferograms.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://olejorik.github.io/PhaseFromInterferograms.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/olejorik/PhaseFromInterferograms.jl",
)
