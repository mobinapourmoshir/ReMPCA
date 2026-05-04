# Plot Method for fdClass Objects

Generates a line plot of the functional data stored in an object of
class `fdClass`.

## Usage

``` r
# S3 method for class 'fdClass'
plot(x, ...)
```

## Arguments

- x:

  An object of class `fdClass`.

- ...:

  Additional graphical parameters passed to plotting functions.

## Value

No return value. This function is called for its side effect (plot).

## Details

This function uses [`matplot`](https://rdrr.io/r/graphics/matplot.html)
to visualize each observation (row) as a separate curve. It provides a
quick overview of the functional data structure stored in the `fdClass`
object.

## Examples

``` r
fd_obj <- fdClass(matrix(sin(1:100 / 10), nrow = 10, ncol = 10))
plot(fd_obj)
#> Warning: There should be an equal number of grid points and columns. 'argval' is set to NULL!
#> Warning: There should be an equal number of grid points and columns. 'argval' is set to NULL!
#> Warning: There should be an equal number of grid points and columns. 'argval' is set to NULL!
#> Warning: There should be an equal number of grid points and columns. 'argval' is set to NULL!
#> Warning: There should be an equal number of grid points and columns. 'argval' is set to NULL!
#> Warning: There should be an equal number of grid points and columns. 'argval' is set to NULL!
#> Warning: There should be an equal number of grid points and columns. 'argval' is set to NULL!
#> Warning: There should be an equal number of grid points and columns. 'argval' is set to NULL!
#> Warning: There should be an equal number of grid points and columns. 'argval' is set to NULL!
#> Warning: There should be an equal number of grid points and columns. 'argval' is set to NULL!

```
