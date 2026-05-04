# Coerce an object of class 'rdClass' or 'imgClass' to class 'fdClass'.

Coerce an object of class 'rdClass' or 'imgClass' to class 'fdClass'.

## Usage

``` r
as.fdClass(
  x,
  Sparsity_parameter = NULL,
  Smoothing_parameter = NULL,
  argval = NULL
)
```

## Arguments

- x:

  An object of class 'rdClass' or 'imgClass'.

- Sparsity_parameter:

  Optional sparsity parameter to assign. If `NULL`, defaults are used.

- Smoothing_parameter:

  Optional smoothing parameter to assign. If `NULL`, defaults are used.

- argval:

  Optional vector of grid points for columns. If `NULL`, defaults are
  used.

## Value

An object of class 'fdClass' with user-specified or inherited
regularization parameters.

## Examples

``` r
img_object <- imgClass(
  image = list(matrix(rnorm(100), nrow = 10), matrix(rnorm(100), nrow = 10)),
  argval = NULL,
  Smoothing_parameter = NULL,
  Sparsity_parameter = 0
)
newfd <- as.fdClass(img_object)
newfd
#> Functional Data (fdClass) Object
#> -----------------------------------
#> Dimensions: 2 x 100
#> Smoothing Parameter: 9.313226e-10 1.379661e-08 2.043829e-07 ... 2.160119 32 
#> Sparsity Parameter: 0 
#> GridPoints_v: 0.01 0.02 0.03 ... 0.99 1 
#> -----------------------------------
#> First few rows and columns of the data:
#>              V1         V2        V3         V4          V5
#> [1,] -1.4000435 -0.5536994 0.4681544  0.9353632  0.07003485
#> [2,] -0.3872136 -0.2088827 2.0393693 -1.1461997 -0.16267634
```
