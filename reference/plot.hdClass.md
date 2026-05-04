# Plot Method for hdClass Objects

Visualizes each variable in a hybrid data object (`hdClass`) using an
appropriate plot style:

- Functional data are shown using lines (`matplot(..., type = "l")`).

- Regular data are shown using solid dots
  (`matplot(..., type = "p", pch = 16)`).

- Image data are plotted using
  [`image()`](https://rdrr.io/r/graphics/image.html) if they originated
  from matrices.

## Usage

``` r
# S3 method for class 'hdClass'
plot(x, ...)
```

## Arguments

- x:

  An object of class `hdClass`.

- ...:

  Additional graphical parameters passed to the plotting functions.

## Value

No return value. Called for its side effect (producing plots).

## Examples

``` r
fd_obj <- fdClass(matrix(rnorm(100), 10, 10))
rd_obj <- rdClass(matrix(rnorm(100), 10, 10))
hd_obj <- hdClass(list(fd_obj, rd_obj))
plot(hd_obj)

```
