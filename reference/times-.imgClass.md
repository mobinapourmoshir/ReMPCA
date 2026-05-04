# Multiply a `imgClass` Object by a Scalar

Multiplies each matrix (image) in an `imgClass` object by a scalar. All
attributes and class structure are retained.

## Usage

``` r
# S3 method for class 'imgClass'
e1 * e2
```

## Arguments

- e1:

  A scalar numeric value or an `imgClass` object.

- e2:

  An `imgClass` object or a scalar numeric value.

## Value

A new `imgClass` object with each image scaled by the scalar value.

## Examples

``` r
img <- imgClass(image = list(matrix(1:9, 3, 3)))
img_scaled <- 2 * img
```
