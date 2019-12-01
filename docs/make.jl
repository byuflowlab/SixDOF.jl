using Documenter, SixDOF

makedocs(;
    modules=[SixDOF],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
    ],
    repo="https://github.com/andrewning/SixDOF.jl/blob/{commit}{path}#L{line}",
    sitename="SixDOF.jl",
    authors="Andrew Ning <aning@byu.edu>",
    # assets=String[],
)
