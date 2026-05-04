# Multiply a `hdClass` Object by a Scalar

Performs element-wise multiplication between a scalar and a `hdClass`
object. All attributes and class information are preserved.

## Usage

``` r
# S3 method for class 'hdClass'
e1 * e2
```

## Arguments

- e1:

  A scalar numeric value or a `hdClass` object.

- e2:

  A `hdClass` object or a scalar numeric value.

## Value

A new `hdClass` object with elements scaled by the scalar value, and
original attributes retained.

## Examples

``` r
obj <- hdClass(list(fdClass(matrix(1:10, ncol = 2))))
2 * obj
#> Hybrid Data (hdClass) Object
#> ===================================
#> Dimensions           : 5 rows x 2 columns
#> Number of Variables  :  1 
#> Columns per Variable :  2 
#> -----------------------------------
#> Grid Points (Rows)   : 0.2, 0.4, 0.6, 0.8, 1 
#> Grid Points (Columns):
#>   - Variable 1: 0.5, 1 
#> Smoothing Parameter (Row): 0 
#> Smoothing Parameters (Col):
#>   - Variable 1: 0 
#> Sparsity Parameter (Row): 0 
#> Sparsity Parameters (Col):
#>   - Variable 1: 0 
#> ===================================
#> First few rows and columns of the data:
#>      V1 V2
#> [1,]  2 12
#> [2,]  4 14
#> [3,]  6 16
#> [4,]  8 18
#> [5,] 10 20
```
