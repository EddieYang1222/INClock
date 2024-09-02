# INClock
INClock (**I**ntrinsic **N**oise based cellular aging **Clock**) provides a robust and scalable measure for cellular aging and senescence by estimating intrinsic transcriptional noise from single-cell RNA sequencing data. By using this package, you are capable of:

1. Compute a manifold that approximates the true state of cells using two methods
  + Decompose cellular variation using [SAVER](https://mohuangx.github.io/SAVER/)
  + Average across k nearest neighbors (more efficient)
2. Estimate the dispersion parameters for each cell types / age groups using [an updated version of SAVER](https://github.com/EddieYang1222/SAVER_D)

## Installation

The most recent version of INClock can be installed from [GitHub](https://github.com/EddieYang1222/INClock) with:

``` r
# install.packages("devtools")
devtools::install_github("EddieYang1222/INClock")
```

R version `> 4.2.0` is recommended.

## Contact

For any questions regarding the package, please contact Yilin Yang (yang1222@wharton.upenn.edu) or Nancy Zhang (nzh@wharton.upenn.edu).

<!-- badges: start -->
<<<<<<< HEAD
[![R-CMD-check](https://github.com/EddieYang1222/INClock/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EddieYang1222/INClock/actions/workflows/R-CMD-check.yaml)
>>>>>>> 4a2580aab71d0b58fb23867b389620f3d1053ff2
<!-- badges: end -->
