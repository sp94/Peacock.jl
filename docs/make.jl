using Documenter, Peacock

makedocs(
    sitename="Peacock.jl Documentation",
    pages = [
        "index.md",
        "Tutorials" => [
            "tutorials/getting_started.md",
        ],
        "How-to guides" => [
            "how-tos/zoo.md",
            "how-tos/wilson_loops.md",
        ],
        "reference/index.md",
        "contributing.md",
    ]
)

deploydocs(
	repo="github.com/sp94/Peacock.jl.git")
