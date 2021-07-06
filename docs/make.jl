using Documenter, DocumenterLaTeX

makedocs(
	sitename="MimiFUND.jl",
	pages=[
		"Home" => "index.md",
		"Science" => "science.md",
		"Tables" => "tables.md"]
)

deploydocs(
    repo="github.com/fund-model/MimiFUND.jl.git"
)
