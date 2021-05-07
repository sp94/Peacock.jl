using Documenter

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
            "how-tos/gpu.md",
        ],
        "reference/index.md",
        "contributing.md",
    ]
)

deploydocs(
    devbranch = "master",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
    repo="github.com/sp94/Peacock.jl.git",
)
