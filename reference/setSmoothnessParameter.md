# Set Smoothness Tuning Parameter

Assigns a smoothing tuning parameter to an object of class `rdClass`,
`fdClass`, `imgClass`, or `hdClass`.

## Usage

``` r
setSmoothnessParameter(obj, Smoothing_parameter)
```

## Arguments

- obj:

  An object of class `rdClass`, `fdClass`, `imgClass`, or `hdClass`.

- Smoothing_parameter:

  A numeric value or vector representing the smoothing parameter(s) to
  assign to the object.

## Value

The modified object with updated `Smoothing_parameter` attribute.

## Details

This function updates the `"Smoothing_parameter"` attribute of the given
object. Since R uses copy-on-modify semantics for S3 objects, users must
reassign the object after calling this function:


      obj <- setSmoothnessParameter(obj, c(0.01, 0.1, 1))

Without reassignment, the original object remains unchanged.

## Examples

``` r
fd_obj <- fdClass(matrix(rnorm(100), 10, 10))
fd_obj <- setSmoothnessParameter(fd_obj, c(0.01, 0.1, 1))
attr(fd_obj, "Smoothing_parameter")
#> [1] 0.01 0.10 1.00
```
