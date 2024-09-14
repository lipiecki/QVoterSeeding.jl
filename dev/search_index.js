var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/#Q-Voter-Model","page":"Functions","title":"Q-Voter Model","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"qvoter","category":"page"},{"location":"functions/#QVoterSeeding.qvoter","page":"Functions","title":"QVoterSeeding.qvoter","text":"qvoter(g::SimpleGraph{<:Integer}, active::BitVector, q::Integer; ε::AbstractFloat = 0.0, flexible::BitVector = trues(nv(g)), replace::Bool = false, warmup::Integer = 10^5, averaging::Integer = 10^3)\n\nRun a Monte Carlo simulation of the q-voter model.\n\nArguments\n\ng: Graph representing the structure of the agent system.\nactive: BitVector representing states of the agents (the vector is mutated during simulation).\nq: Number of neighbors in the influence group.\n\nKeyword arguments\n\nε = 0.0: Probability of flipping state when q-panel is not unanimous.\nflexible = trues(nv(g)): BitVector represeting flexibility of the agents, i.e., whether an agent changes its state during simulation.\nreplace = false: If true, cosntruct the q-panel by drawing neighbors with replacement, otherwise draw without replacement.\nwarmup = 10^5: Number of Monte Carlo steps calculated during the simulation warm up.\naveraging = 10^3: Number of Monte Carlo steps over which the concentration is averaged.\n\nReturn the final concentration of active agents. If the system reaches an absorbing state, the final concentration is returned immediately.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Seeding","page":"Functions","title":"Seeding","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"seeding","category":"page"},{"location":"functions/#QVoterSeeding.seeding","page":"Functions","title":"QVoterSeeding.seeding","text":"seeding(network::Union{Symbol, String}, strategy::Symbol, zealots::Bool, q::Integer, budgets::AbstractVector{<:AbstractFloat}, rng::AbstractRNG = Random.default_rng())\n\nRun network seeding experiments with q-voter dynamics with for a specified seeding strategy and seeding budgets. If zealots = true, seeds will act as zealots.\n\nArguments network and rng are passed to loadnetwork.\n\nAvailable options for strategy:\n\n:hd for high degree\n:pr for PageRank\n:cc for complex centrality\n:onehop for one-hop\n:random for random.\n\nReturn a vector of final concentrations of active vertices for each seeding budget.    \n\n\n\n\n\n","category":"function"},{"location":"functions/#Networks","page":"Functions","title":"Networks","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"loadnetwork\naverage_degree\naverage_clustering\naverage_pathlen\nprune!\naddneighbors!","category":"page"},{"location":"functions/#QVoterSeeding.loadnetwork","page":"Functions","title":"QVoterSeeding.loadnetwork","text":"loadnetwork(networkname::Symbol, rng::AbstractRNG = Random.default_rng(); preprocess::Bool = true)\n\nLoad network data from the networks directory or from the specified filepath.\n\nAvailable options for network:\n\n:snap for SNAP Facebook network\n:politicians for GEMSEC Facebook page network of politicians\n\"...\" for network at filepath \"...\"; the file is assumed to be an edge list in a CSV format with node indexing starting from 0.\n\nOptional argument rng specifies a random number generator used for preprocessing of the network data. Keyword argument preprocess specifies whether to perform preprocessing (prune! and addneighbors!).\n\nReturn a graph representing the selected network.\n\n\n\n\n\n","category":"function"},{"location":"functions/#QVoterSeeding.average_degree","page":"Functions","title":"QVoterSeeding.average_degree","text":"average_degree(network::Union{Symbol, String}, m::Integer, rng::AbstractRNG = Random.default_rng())\n\nCalculate the average node degree of a raw network and an ensemble of m preprocessed networks.\n\nAvailable options for network:\n\n:snap for SNAP Facebook network\n:politicians for GEMSEC Facebook page network of politicians\n\"...\" for network file at path \"...\".\n\nOptional argument rng specifies a random number generator used for generating an ensemble of preprocessed networks.\n\nReturn a tuple of the average node degree for a raw network and for a network ensemble.\n\n\n\n\n\n","category":"function"},{"location":"functions/#QVoterSeeding.average_clustering","page":"Functions","title":"QVoterSeeding.average_clustering","text":"average_clustering(network::Union{Symbol, String}, m::Integer, rng::AbstractRNG = Random.default_rng())\n\nCalculate the average clustering coefficient of a raw network and an ensemble of m preprocessed networks.\n\nAvailable options for network:\n\n:snap for SNAP Facebook network\n:politicians for GEMSEC Facebook page network of politicians\n\"...\" for network file at path \"...\".\n\nOptional argument rng specifies a random number generator used for generating an ensemble of preprocessed networks.\n\nReturn a tuple of the average clustering coefficient for a raw network and for a network ensemble.\n\n\n\n\n\n","category":"function"},{"location":"functions/#QVoterSeeding.average_pathlen","page":"Functions","title":"QVoterSeeding.average_pathlen","text":"average_pathlen(network::Union{Symbol, String}, m::Integer, rng::AbstractRNG = Random.default_rng())\n\nCalculate the average shortest path length of a raw network and an ensemble of m preprocessed networks.\n\nAvailable options for network:\n\n:snap for SNAP Facebook network\n:politicians for GEMSEC Facebook page network of politicians\n\"...\" for network file at path \"...\".\n\nOptional argument rng specifies a random number generator used for generating an ensemble of preprocessed networks.\n\nReturn a tuple of the average shortest path length for a raw network and for a network ensemble.\n\n\n\n\n\n","category":"function"},{"location":"functions/#QVoterSeeding.prune!","page":"Functions","title":"QVoterSeeding.prune!","text":"prune!(g::SimpleGraph{<:Integer})\n\nModify graph g by removing self-loops and isolated vertices.\n\n\n\n\n\n","category":"function"},{"location":"functions/#QVoterSeeding.addneighbors!","page":"Functions","title":"QVoterSeeding.addneighbors!","text":"addneighbors!(g::SimpleGraph{<:Integer}, rng::AbstractRNG = Random.default_rng())\n\nModify graph g by adding a new neighbor to every single-neighbor-vertex via triad formation.\n\nOptional argument rng specifies a random number generator.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Complex-Centrality","page":"Functions","title":"Complex Centrality","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"complexcentrality\ncompetitive_complexcentrality","category":"page"},{"location":"functions/#QVoterSeeding.complexcentrality","page":"Functions","title":"QVoterSeeding.complexcentrality","text":"complexcentrality(g::SimpleGraph{<:Integer}, T::Vector{<:Integer}[, v::Integer])\n\nCalculate the complex centrality for graph g with threshold values specifified by vector T.\n\nReturn a vector of complex centralities calculated for each vertex. If v is specified, calculate and return only the complex centrality of vertex v.\n\n\n\n\n\n","category":"function"},{"location":"functions/#QVoterSeeding.competitive_complexcentrality","page":"Functions","title":"QVoterSeeding.competitive_complexcentrality","text":"competitive_complexcentrality(g::SimpleGraph{<:Integer}, q::Integer[, v::Integer])\n\nCalculate the complex centrality for graph g with threshold values adjusted to q-voter dynamics, according to: T_i = maxq lceil k_i2rceil, where k_i is the degree of vertex i.\n\nReturn a vector of complex centralities calculated for each vertex. If v is specified, calculate and return only the complex centrality of vertex v.\n\n\n\n\n\n","category":"function"},{"location":"functions/#One-hop","page":"Functions","title":"One-hop","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"onehop","category":"page"},{"location":"functions/#QVoterSeeding.onehop","page":"Functions","title":"QVoterSeeding.onehop","text":"onehop(g::SimpleGraph{<:Integer}, seedsize::Integer)\n\nConstruct a seed set of size seedsize on graph g using one-hop strategy.\n\nReturn a vector of seed vertices.\n\n\n\n\n\n","category":"function"},{"location":"#QVoterSeeding.jl","page":"Home","title":"QVoterSeeding.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Julia package for network seeding simulations under q-voter dynamics","category":"page"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is supplementary to the paper Strategic Promotional Campaigns for Sustainable Behaviors: Maximizing Influence in Competitive Complex Contagions, if you find it useful in your research, please consider adding the following citation:","category":"page"},{"location":"","page":"Home","title":"Home","text":"@inproceedings{lipiecki2024,\n    author = {Lipiecki, Arkadiusz},\n    year = {2024},\n    title = {Strategic Promotional Campaigns for Sustainable Behaviors: Maximizing Influence in Competitive Complex Contagions},\n    editor = {Franco, Leonardo\n              and de Mulatier, Cl{\\'e}lia\n              and Paszynski, Maciej\n              and Krzhizhanovskaya, Valeria V.\n              and Dongarra, Jack J.\n              and Sloot, Peter M. A.},\n    booktitle = {Computational Science -- ICCS 2024},\n    publisher = {Springer Nature Switzerland},\n    address = {Cham},\n    pages = {62--70},\n    isbn = {978-3-031-63759-9},\n    doi = {https://doi.org/10.1007/978-3-031-63759-9_8}\n}","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the package, use the following lines in Julia REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg\njulia> Pkg.add(url = \"https://github.com/lipiecki/QVoterSeeding.jl\")","category":"page"},{"location":"#Datasets","page":"Home","title":"Datasets","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The following network datasets are provided with the package:","category":"page"},{"location":"","page":"Home","title":"Home","text":"SNAP Facebook Network\nGEMSEC Facebook Page Network of Politicians","category":"page"},{"location":"#Example","page":"Home","title":"Example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This code snippet performs network seeding experiments on the SNAP Facebook network using the one-hop seeding strategy within the 2-voter model with zealot (inflexible) seeds:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using QVoterSeeding\n\nnetwork = :snap\nstrategy = :onehop\nzealots = true\nq = 2\nbudgets = 0.025:0.025:0.5\n\nresults = seeding(network, strategy, zealots, q, budgets)","category":"page"},{"location":"","page":"Home","title":"Home","text":"results is a vector of final concentrations of active nodes obtained for each seeding budget","category":"page"}]
}
