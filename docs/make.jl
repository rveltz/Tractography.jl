ENV["GKSwstype"] = "100" 

using Documenter, CairoMakie, Tractography

# to display progress
ENV["JULIA_DEBUG"] = Documenter

makedocs(
	modules = [Tractography, isdefined(Base, :get_extension) ?
    Base.get_extension(Tractography, :MakieExt) :
    Tractography.MakieExt],
	doctest = false,
	sitename = "Tractography in Julia",
	format = Documenter.HTML(edit_link = "master", collapselevel = 1),
	# format = Documenter.LaTeX(),
	warnonly = true,
	draft = false,
	pagesonly = true, # do not compile what is not in pages =
	authors = "Samuel Deslauriers-Gauthier and Romain Veltz",
	pages = Any[
		"Home" => "index.md",
		"Getting started" => "getstart.md",
		"Tutorials" => [
				# "basic.md",
				"fcup.md",
				"gpu.md",
		],
		"Algorithms" => "algos.md",
		"Plotting" => "plot.md",
		"Library" => "library.md"
	],
	remotes  = nothing,
	)