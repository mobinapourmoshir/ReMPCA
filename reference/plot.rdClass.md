# Plot Method for rdClass Objects

Produces a scatter plot of regular data stored in an `rdClass` object.
Each column is plotted as a sequence of solid points.

## Usage

``` r
# S3 method for class 'rdClass'
plot(x, ...)
```

## Arguments

- x:

  An object of class `rdClass`.

- ...:

  Additional graphical parameters passed to plotting functions.

## Value

No return value. This function is called for its side effect (plot).

## Details

This method visualizes the regular (non-functional) data in the
`rdClass` object. It shows the values in each column as solid dots,
which is useful for examining patterns across observations or variables.

## Examples

``` r
rd_obj <- rdClass(matrix(rnorm(100), nrow = 10, ncol = 10))
plot(rd_obj)

```
