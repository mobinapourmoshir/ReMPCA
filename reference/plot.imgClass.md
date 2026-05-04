# Plot Method for imgClass Objects

Visualizes a list of image matrices stored in an `imgClass` object. Each
image is shown one at a time. The user is prompted to press Enter or
click on the graphics window to advance to the next image.

## Usage

``` r
# S3 method for class 'imgClass'
plot(x, ...)
```

## Arguments

- x:

  An object of class `imgClass`, created using
  [`imgClass()`](https://mobinapourmoshir.github.io/ReMPCA/reference/imgClass.md).

- ...:

  Additional graphical parameters passed to plotting functions.

## Value

No return value. Called for its side effect (plotting).

## Examples

``` r
img_list <- list(matrix(rnorm(100), 10, 10),
                 matrix(runif(100), 10, 10))
img_obj <- imgClass(img_list)
plot(img_obj)

#> Press [Enter] or click to continue...

```
