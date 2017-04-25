# mutools3D 

The **mutools3D** package implements R functions for mass univariate analysis of three-dimensional phenotypes. 

## Installation

Currently this package is only available on GitHub and requires R (>= 3.3.2). For installation, the easiest way is using `devtools`. If you don’t have `devtools` installed, please start by typing:

```r
install.packages("devtools")
library(devtools)
```

Then, to install the package just use `install_github`function.

```r
install_github("UK-Digital-Heart-Project/mutools3D", build_vignettes = TRUE)
```

## Documentation
An introductory vignette is available from within R (and also [here](inst/doc/gettingStarted.pdf)). Load the package and then use the vignette function:

```r
library(mutools3D)
vignette("gettingStarted", package = "mutools3D")
```

The code has been fully documented and function descriptions are available from within R, i.e. `help(mur)`.


