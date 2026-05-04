# Hybrid Data Constructor (`hdClass`)

Constructs an object of class `hdClass`, representing hybrid data
composed of multiple functional and/or raw variables.

## Usage

``` r
hdClass(hdlist, argval = NULL, Smoothing_parameter = 0, Sparsity_parameter = 0)
```

## Arguments

- hdlist:

  A list of `fdClass`, `rdClass`, or `imgClass` objects.

- argval:

  Optional numeric vector of grid points along the rows. If `NULL`, it
  defaults to a uniform grid over 0 and 1.

- Smoothing_parameter:

  Smoothing parameter(s) for the rows:

  - If 0, no smoothing is applied.

  - If a numeric value or vector, it is used directly.

  - If `NULL`, defaults to `2^seq(-30, 5, length.out = 10)` for tuning.

- Sparsity_parameter:

  Sparsity parameter(s) for the rows:

  - If 0, no sparsity is applied.

  - If a numeric vector, it is used for tuning.

  - If `NULL`, a suitable default sequence is generated.

## Value

An object of class `hdClass` (a matrix) with several hybrid-aware
attributes.

## Details

The `hdClass` is an S3 class for representing *hybrid data*, consisting
of a list of functional (`fdClass`), raw (`rdClass`), or image
(`imgClass`) data objects.

- The input list (`hdlist`) can contain a mix of `fdClass`, `rdClass`,
  and `imgClass` objects.

- All objects must have the **same number of observations** (i.e., same
  number of rows).

- Image objects (`imgClass`) are interpreted as either functional or raw
  based on their internal class inheritance.

The resulting object is a matrix with column-wise concatenation of the
input matrices and includes the following attributes:

- `GridPoints_u`: Row grid points.

- `Smoothing_parameter`: Row smoothing parameter(s).

- `Sparsity_parameter`: Row sparsity parameter(s).

- `n_var`: Number of variables (i.e., number of objects in `hdlist`).

- `ncol`: Number of columns in each variable.

- `Smoothing_parameter_col`: List of smoothing parameters for each
  variable.

- `Sparsity_parameter_col`: List of sparsity parameters for each
  variable.

- `GridPoints_v`: List of grid points (columns) for each variable.

- `variable_types`: Character vector of type labels for each variable:
  `"hd"` (functional) or `"rd"` (raw).

## Examples

``` r
fd_obj <- fdClass(matrix(rnorm(100), nrow = 10))
rd_obj <- rdClass(matrix(rnorm(100), nrow = 10))

hd_obj <- hdClass(
  hdlist = list(fd_obj, rd_obj),
  argval = seq(0, 1, length.out = 10),
  Smoothing_parameter = 0.5,
  Sparsity_parameter = 2
)

print(hd_obj)
#> Hybrid Data (hdClass) Object
#> ===================================
#> Dimensions           : 10 rows x 20 columns
#> Number of Variables  :  2 
#> Columns per Variable :  10 | 10 
#> -----------------------------------
#> Grid Points (Rows)   : 0, 0.111111111111111, 0.222222222222222 , ...,  0.888888888888889, 1 
#> Grid Points (Columns):
#>   - Variable 1: 0.1, 0.2, 0.3 , ...,  0.9, 1 
#>   - Variable 2: NULL
#> Smoothing Parameter (Row): 0.5 
#> Smoothing Parameters (Col):
#>   - Variable 1: 0 
#>   - Variable 2: 0 
#> Sparsity Parameter (Row): 2 
#> Sparsity Parameters (Col):
#>   - Variable 1: 0 
#>   - Variable 2: 0 
#> ===================================
#> First few rows and columns of the data:
#>                V1         V2          V3         V4         V5
#> [1,]  0.893165020  2.2729548  0.88986536 -0.1865538 -0.5296437
#> [2,] -0.376555496  0.1360582 -1.48281332 -0.8122754 -0.6812732
#> [3,]  0.605884808 -1.9990391  0.44575035 -1.6405817 -0.2024476
#> [4,] -0.004874726 -0.4205009  1.36977586  0.5079225  1.6844957
#> [5,] -0.520796373 -0.3784074 -0.02011003  1.7543370 -1.0337732
attr(hd_obj, "variable_types")  # Shows "hd", "rd", etc.
#> [1] "hd" "rd"
```
