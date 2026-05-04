# Set Sparsity Tuning Parameter

Assigns a sparsity tuning parameter to an object of class `rdClass`,
`fdClass`, `imgClass`, or `hdClass`.

## Usage

``` r
setSparsityParameter(obj, Sparsity_parameter)
```

## Arguments

- obj:

  An object of class `rdClass`, `fdClass`, `imgClass`, or `hdClass`.

- Sparsity_parameter:

  A numeric vector of non-negative integers indicating sparsity tuning
  levels. If `NULL`, a default sequence will be generated based on the
  object's dimensions.

## Value

The modified object with updated sparsity parameters.

## Details

This function returns a modified version of the input object with an
updated `"Sparsity_parameter"` attribute. Due to R's copy-on-modify
semantics for S3 objects, the user must reassign the object:


      obj <- setSparsityParameter(obj, c(0, 2, 4))

The object will not be updated in-place unless reassigned.

## Examples

``` r
fd_obj <- fdClass(matrix(rnorm(100), 10, 10))
fd_obj <- setSparsityParameter(fd_obj, 0:5)
attr(fd_obj, "Sparsity_parameter")
#> [1] 0 1 2 3 4 5
```
