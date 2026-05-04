# Coerce an object of class 'rdClass', 'fdClass', or 'imgClass' to class 'hdClass'

This function wraps a single regular or functional object into a hybrid
data object.

## Usage

``` r
as.hdClass(
  x,
  Smoothing_parameter = NULL,
  Sparsity_parameter = NULL,
  argval = NULL
)
```

## Arguments

- x:

  An object of class `'rdClass'`, `'fdClass'`, or `'imgClass'`.

- Smoothing_parameter:

  Optional smoothing parameter for the row direction. If `NULL`, taken
  from input attributes or defaulted in
  [`hdClass()`](https://mobinapourmoshir.github.io/ReMPCA/reference/hdClass.md).

- Sparsity_parameter:

  Optional sparsity parameter for the row direction. If `NULL`, taken
  from input attributes or defaulted in
  [`hdClass()`](https://mobinapourmoshir.github.io/ReMPCA/reference/hdClass.md).

- argval:

  Optional grid points for rows. If `NULL`, taken from input attributes
  or defaulted.

## Value

An object of class `'hdClass'`.

## Details

This function returns a modified version of the input object as an
`hdClass` object. Due to R's copy-on-modify semantics for S3 objects,
users must reassign the result to retain changes. For example:
`x <- as.hdClass(x)`.

## Examples

``` r
fd_obj <- fdClass(matrix(rnorm(100), 10, 10))
hd_obj <- as.hdClass(fd_obj, Sparsity_parameter = 0:5)
attr(hd_obj, "Sparsity_parameter")
#> [1] 0 1 2 3 4 5
```
