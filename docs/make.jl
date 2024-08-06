using Documenter, DocumenterVitepress

using Subzero

makedocs(;
    modules=[Subzero],
    authors="Skylar Gering",
    repo="https://github.com/Caltech-OCTO/Subzero.jl",
    sitename="Subzero.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/Caltech-OCTO/Subzero.jl",
        devurl = "dev",
        deploy_url = "Caltech-OCTO.github.io/Subzero.jl",
    ),
    pages=[
        "Home" => "index.md",
    ],
    warnonly = true,
)

deploydocs(;
    repo="https://github.com/Caltech-OCTO/Subzero.jl",
    push_preview=true,
)
