using Documenter, QVoterSeeding

makedocs(
    sitename = "QVoterSeeding.jl",
    format = Documenter.HTML(
        repolink = "https://github.com/lipiecki/QVoterSeeding.jl",
        inventory_version = ""),
    pages = [
        "Home" => "index.md",
        "Functions" => "functions.md"])
deploydocs(
    repo = "github.com/lipiecki/QVoterSeeding.jl.git",
    branch = "gh-pages")
