
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sparsemixwishart

The goal of sparsemixwishart is to provide a framework for sparse
model-based clustering of covariance objects. The approach assumes that
the sample of covariance matrices are iid from a mixture of Wishart
distribution and the parameters have different degree of sparsity across
clusters, allowing to induce parsimony in a flexible manner. Estimation
of the model relies on the maximization of a penalized likelihood, with
a specifically tailored covariance graphical lasso penalty.

This repository is associated with the paper Cappozzo, Casa (2024+)
*Model-based clustering for covariance matrices via penalized Wishart
mixture models* <!-- [FIXME arxiv link](FIXME_arxiv_link) -->

## Installation

You can install the development version of sparsemixwishart from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AndreaCappozzo/sparsemixwishart")
```
