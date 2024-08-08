using Documenter, Literate

using Subzero

# Converting tutorial to markdown
Literate.markdown(
    joinpath(@__DIR__, "src", "tutorial.jl"), joinpath(@__DIR__, "src");
    credit = false
)

format = Documenter.HTML(;
    repolink = "https://github.com/Caltech-OCTO/Subzero.jl",
    canonical = "https://Caltech-OCTO.github.io/SubzeroDocumentation/stable",
    mathengine = MathJax3(),
    size_threshold = 5*10^5,  # 500 KiB
)

makedocs(;
    modules=[Subzero],
    authors="Skylar Gering and contributers",
    sitename="Subzero.jl",
    format,
    pages=[
        "Introduction" => "index.md",
        "Tutorial" => "tutorial.md",
        "API Reference" => "api.md",
        "Developer Documentation" => "devdocs.md"
    ],
    warnonly = true,
)

deploydocs(;
    repo="https://github.com/Caltech-OCTO/Subzero.jl",
)
