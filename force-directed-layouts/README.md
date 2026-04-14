**This folder generates force-directed layouts for the Omi32 combinatorial library in Tharp et al, 2026**
- The force-directed layout sets each evolutionary intermediate in the library (each unique genotype) as a node, and each node is connected by edges to each of its single mutation neighbors
- Edges are set as springs, which have 'weights' (i.e., spring constants) set such that two connected nodes pull on each other proportional to the change in phenotype between them (i.e., connected nodes with small changes in phenotype tend to cluster together)
- this version of the browser uses the BA1 affinities (phenotype) to generate a force-directed layout and a pyarrow file (.pyarrow) which are used to construct the interactive web browser, but any of the phenotypes (or fitness values from `pathway_inference`) can be used to make these graphs and display in your own interactive web browser!
- see https://amphilli.github.io/Omi32_browser_git/ for the published interactive web-browser

** notes:**
- running this notebook requires the igraph package
- for simplicity, we've included `environment.yml` in this folder for running the notebook
- set this folder as your current working directory, then enter `conda env create -f environment.yml` and use `conda activate force-directed-layout` (alternatively, use your favorite code editor and select the `structural-analysis` conda environment that you just installed)