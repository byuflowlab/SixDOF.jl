using Documenter, SixDOF

makedocs(;
    modules=[SixDOF],
    format=Documenter.HTML(),
    pages=[
        "Guide" => "guide.md",
        "Theory" => "theory.md",
        "Reference frames" => "frames.md",
    ],
    sitename="SixDOF.jl",
    authors="Andrew Ning <aning@byu.edu>",
    checkdocs=:exports,
    warnonly=[:missing_docs],
)

deploydocs(
    repo = "github.com/byuflowlab/SixDOF.jl.git",
)