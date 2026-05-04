# Coerce an object of class 'fdClass' or 'imgClass' to class 'rdClass'

This function strips the smoothing and grid-related structure from a
functional object and converts it to a regular data object of class
`rdClass`.

## Usage

``` r
as.rdClass(x, Sparsity_parameter = NULL)
```

## Arguments

- x:

  An object of class `'fdClass'` or `'imgClass'`.

- Sparsity_parameter:

  Optional sparsity parameter to override the original. If `NULL`, the
  existing sparsity parameter (if any) is inherited.

## Value

An object of class `rdClass` with smoothing and grid attributes removed,
and sparsity parameter preserved or overridden.

## Examples

``` r
img_object <- imgClass(image = list(matrix(rnorm(100), nr = 50),
                                    matrix(rnorm(100), nr = 50)),
                       argval = NULL,
                       Smoothing_parameter = NULL,
                       Sparsity_parameter = 0)

newrd <- as.rdClass(img_object, Sparsity_parameter = 1:10)
attr(newrd, "Sparsity_parameter")
#>  [1]  1  2  3  4  5  6  7  8  9 10
```
