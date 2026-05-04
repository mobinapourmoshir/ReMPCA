# Element-wise Subtraction for Hybrid Data Objects

Performs element-wise subtraction of two objects of class `hdClass`,
`fdClass`, `rdClass`, or `imgClass`, assuming they have identical
dimensions. This operator is useful for iterative model fitting,
residual computation, or gradient-based updates in hybrid data
decomposition.

## Usage

``` r
# S3 method for class 'hd'
obj1 - NULL
```

## Arguments

- obj1:

  An object of class `hdClass`, `fdClass`, `rdClass`, or `imgClass`.

- obj2:

  Another object of the same class as `obj1`. If `NULL`, the function
  returns `obj1`.

## Value

An object of the same class as `obj1` and `obj2`, representing the
element-wise difference.

## Details

The dimensions of the two input objects must match exactly. Attributes
from `obj1` are preserved.

## Examples

``` r
fd1 <- fdClass(matrix(1:9, 3, 3))
fd2 <- fdClass(matrix(1, 3, 3))
fd_diff <- fd1 - fd2
print(fd_diff)
#> Functional Data (fdClass) Object
#> -----------------------------------
#> Dimensions: 3 x 3
#> Smoothing Parameter: 0 
#> Sparsity Parameter: 0 
#> GridPoints_v: 0.3333333 0.6666667 1 
#> -----------------------------------
#> First few rows and columns of the data:
#>      V1 V2 V3
#> [1,]  0  3  6
#> [2,]  1  4  7
#> [3,]  2  5  8
```
