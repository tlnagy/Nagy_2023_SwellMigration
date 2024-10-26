# Overview

!!! note

    This mini-site is meant to accompany the paper:

    T. L. Nagy, E. Strickland, O. D. Weiner, "Neutrophils actively swell to 
    potentiate rapid migration" *Elife* **2023**, 
    [DOI 10.7554/elife.90551.](https://elifesciences.org/reviewed-preprints/90551)

The goal is facilitate reproducible and reusable research by publishing the code
and data such that the analysis is built using continuous integration to make
sure that the code is valid, functional, and robust to future updates. The
[Figures](@ref) page consists of separate `Literate.jl` scripts used to make the
figures in the paper.

The easiest way to get all the code working on your computer is to:

1. [Install a recent version of Julia](https://julialang.org/downloads/)

2. [Clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
   the [project's Github repo](https://github.com/tlnagy/Nagy_2023_SwellMigration) 
   to your computer

3. Run the following two lines in your terminal:

```sh
julia --project="." -e 'using Pkg; Pkg.instantiate()' 

julia --project="." site/make.jl
```

This will download and install all dependencies and then build this site. You
can then run the following to serve the website:

```sh
julia --project=build -e 'using LiveServer; serve(dir = "site/build")'
```

The server will print out the URL that you can then enter into your browser's
navigation bar to see the website rendered nicely.

!!! tip
    The scripts can be downloaded individually using the 
    ```@raw html
    <img src="https://img.shields.io/badge/download-julia-brightgreen.svg"
    alt="Source code"/>
    ``` or 
    ```@raw html
    <img src="https://img.shields.io/badge/show-nbviewer-579ACA.svg"
    alt="notebook"/> 
    ``` buttons on the top of each page. However, you'll need to make sure that the Julia environment is set up properly and that the data is
    located in the correct relative path compared to the script.

    The data can be downloaded independently:

    > Nagy, Tamas; Strickland, Jack; Weiner, Orion (2023). Data from: Neutrophils actively swell to potentiate rapid migration [Dataset]. Dryad. <https://doi.org/10.7272/Q6NS0S5N>

    and should be placed in a folder called `data/` in the parent folder containing the scripts

This site would not be possible without the fantastic
[`Documenter.jl`](https://documenter.juliadocs.org/stable/), 
[`Literate.jl`](https://github.com/fredrikekre/Literate.jl), and
[`DemoCards.jl`](https://democards.juliadocs.org/stable/) libraries. 

~Tamas