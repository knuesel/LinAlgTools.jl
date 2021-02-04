using LinAlgTools
using Documenter

makedocs(;
    modules=[LinAlgTools],
    authors="Jeremie Knuesel <knuesel@gmail.com> and contributors",
    repo="https://github.com/knuesel/LinAlgTools.jl/blob/{commit}{path}#L{line}",
    sitename="LinAlgTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://knuesel.github.io/LinAlgTools.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/knuesel/LinAlgTools.jl",
)
