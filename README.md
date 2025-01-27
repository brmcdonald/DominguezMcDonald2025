# Enrichment of microbial DNA in plasma to improve pathogen detection in sepsis: a pilot study
Data and analysis code for Dominguez, McDonald et al, Clinical Chemistry 2025

To run analysis:
Download Julia (developed using Julia 1.10.4)   
https://julialang.org/downloads/

Add the julia environment for the analysis in julia  
```julia -e "using Pkg; Pkg.develop(path=\"./env_DominguezMcDonald2025\")"```

Run the primary analysis script   
```julia ./anaylsisCode_DominguezMcDonald2025.jl```

To change the input/output directory paths, edit the paths at the top of the anaylsisCode_DominguezMcDonald2025.jl file.
