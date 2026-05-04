# Multiply a `fdClass` Object by a Scalar

Performs element-wise multiplication between a scalar and a `fdClass`
object. All functional data attributes are preserved in the result.

## Usage

``` r
# S3 method for class 'fdClass'
e1 * e2
```

## Arguments

- e1:

  A scalar numeric value or a `fdClass` object.

- e2:

  A `fdClass` object or a scalar numeric value.

## Value

A new `fdClass` object with elements scaled by the scalar value.

## Examples

``` r
fd <- fdClass(matrix(1:20, 10, 2), Smoothing_parameter = 0.1)
scaled_fd <- 3 * fd
```
