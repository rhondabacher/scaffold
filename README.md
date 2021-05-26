# Scaffold: data-generation simulaton of single-cell RNA-seq data

`Scaffold` is an R package for generating scRNA-seq data by
statistically modeling each step of the experimental data generation
process. Scaffold is able to simulate non-UMI, UMI, and 10X data. 

# Installation

Scaffold can be installed from Github:

```{r}
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")
devtools::install_github("rhondabacher/Scaffold")
```

# Vignette

Details and examples on how to use R package are located in the vignette:

* [Vignette](XXX)


# Citation

A preprint of the Scaffold manuscript is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.10.05.326553v1.abstract).


# Contact/Maintainer

* [Rhonda Bacher](https://www.rhondabacher.com) (rbacher@ufl.edu)  
Department of Biostatistics, University of Florida
