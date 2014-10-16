# Fast SCCA (fscca)

**This project is a work in progress**

This is a port of NIPALS SCCA package written by Lee et. al.
\[[1](http://fafner.meb.ki.se/personal/yudpaw/?page_id=13),
[2](http://www.degruyter.com/view/j/sagmb.2011.10.issue-1/sagmb.2011.10.1.1638/sagmb.2011.10.1.1638.xml)\]
(scca). It is written for speedup in Rcpp.

Currently, only the lasso penalty function is implemented. Also, the only
cross-validation method implemented is the sequential one-dimensional search
described in the paper.

# Installation

Make sure you have the `devtools` package installed as well as `Rcpp` and
`RcppArmadillo`:

```R
install.packages(c("devtools", "Rcpp", "RcppArmadillo"))
library(devtools)
```

You can install the current development version from github using `devtools`:

```R
install_github("pimentel/fscca")
```

The build sets the flag `PKG_CXXFLAGS = "-std=c++11"` in `src/Makevars`. This
is standard for `g++` and `clang`. If this fails on your computer, let me know
(my guess it will likely fail on windows).

The branch `master` should always be stable. Please file a bug report if it is
unstable.

Assuming all goes well, load it like a usual package:

```R
library(fscca)
```

# Usage

The main `R` entry point is the function `fscca()`. It is reasonably well
documented in `R` which you can access by typing `?fscca` or `help('fscca')`.
This function will automatically perform cross-validation as well as compute
multiple components (default is just the first).

If you don't want to perform cross-validation, you can compute the first
canonical vectors by using `nipals_sparse()` which is also documented in `R`.

# TODOs

* Implement other penalty functions (soft-thresholding, SCAD, HL)
* Three matrix extension ( `fscca(X, Y, Z)` )
* Implement full grid search

# Bugs

Please report any bugs to Github:
[http://github.com/pimentel/fast-scca/issues](http://github.com/pimentel/fast-scca/issues)
