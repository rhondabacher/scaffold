https://img.shields.io/github/downloads/rhondabacher/scaffold/total

# Scaffold: data-generation simulation of single-cell RNA-seq data

`Scaffold` is an R package for generating scRNA-seq data by
statistically modeling each step of the experimental data generation
process. Scaffold is able to simulate non-UMI, UMI, and 10X data. 

# Installation

Scaffold can be installed from Github:

```{r}
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")
devtools::install_github("rhondabacher/scaffold")
```

# Vignette

Details and examples on how to use the Scaffold R package are located:

* [Scaffold Vignette](https://www.rhondabacher.com/scaffold-vignette.pdf)


# Citation

Additional details on the Scaffold simulation scheme are available in the manuscript: 

[Rhonda Bacher, Li-Fang Chu, Cara Argus, Jennifer M Bolin, Parker Knight, James A Thomson, Ron Stewart, Christina Kendziorski, Enhancing biological signals and detection rates in single-cell RNA-seq experiments with cDNA library equalization, Nucleic Acids Research, Volume 50, Issue 2, 25 January 2022, Page e12, https://doi.org/10.1093/nar/gkab1071](https://academic.oup.com/nar/article/50/2/e12/6438026).


# Contact/Maintainer

* [Rhonda Bacher](https://www.rhondabacher.com) (rbacher@ufl.edu)  
Department of Biostatistics, University of Florida
