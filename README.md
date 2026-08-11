# GenomicBreedingPlots

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GenomicBreeding.github.io/GenomicBreedingPlots.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GenomicBreeding.github.io/GenomicBreedingPlots.jl/dev/)
[![Build Status](https://github.com/GenomicBreeding/GenomicBreedingPlots.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GenomicBreeding/GenomicBreedingPlots.jl/actions/workflows/CI.yml?query=branch%3Amain)

Non-interactive plotting library for GenomicBreeding.jl

## Dev stuff:

### REPL prelude

```shell
julia --project=. -e 'using Pkg; Pkg.instantiate()' # For a fresh Julia installation
julia --project=. --threads=2,1 --load test/interactive_prelude.jl
```

### Format and test

```shell
time julia --project=. --threads=2 -e "using Pkg; Pkg.update()"
time julia --project=. --threads=2  test/cli_tester.jl
```
