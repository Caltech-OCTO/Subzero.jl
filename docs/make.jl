using Documenter

using Subzero

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
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "API Reference" => "api.md",
    ],
    warnonly = true,
)

deploydocs(;
    repo="https://github.com/Caltech-OCTO/Subzero.jl",
)
