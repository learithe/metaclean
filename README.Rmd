<!-- README.md is generated from README.Rmd. Please edit that file -->
metaclean
==========

Functions for cleaning up metabarcoding datasets


### Motivation

This package contains a set of functions I use to work with metabarcoding data output from [DADA2](https://benjjneb.github.io/dada2/), to remove contamination, convert dada2 sequences into easier-to-work-with sequence IDs, and generate biologist-friendly output tables.

### Expected Input
[Phyloseq](https://joey711.github.io/phyloseq/) objects


### Installation

You can install metaclean from github with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("learithe/metaclean")
```


#### Dependencies

metaclean requires the packages **phyloseq** and **dplyr**. It has only been tested with phyloseq version ‘1.22.3’.
