using Documenter, SixDOF

makedocs(;
    modules=[SixDOF],
    format=Documenter.HTML(),
    pages=[
        "Guide" => "guide.md",
        "Theory" => "theory.md",
    ],
    repo="https://github.com/byuflowlab/SixDOF.jl/blob/{commit}{path}#L{line}",
    sitename="SixDOF.jl",
    authors="Andrew Ning <aning@byu.edu>",
    # assets=String[],
)

deploydocs(
    repo = "github.com/byuflowlab/SixDOF.jl.git",
)