# Hierarchical_WLC
This repository contains code used for generating force-extension curves of the Hierarchical worm-like chain (HWLC) model in the paper "Nonlinear mechanics of human mitotic chromosomes".

The code runs using Julia 1.6 or 1.7, which can be installed at https://julialang.org/downloads/.

Required packages:
  - Roots
  - Plots
  - Distributions
  - MAT

### Summary of files and folders:
1. WLC.jl : calculate extensions for single WLC
2. HWLC.jl : module for calculating extensions of HWLCs
3. example_FOO.jl : generate plots for FOO distribution
4. Experimental_data/ : .mat files with experimental data
5. Plots/ : default path for saving generated plots
