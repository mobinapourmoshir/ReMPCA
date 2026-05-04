# Multiply a `rdClass` Object by a Scalar

Performs element-wise multiplication between a scalar and a `rdClass`
object. Attributes like sparsity settings are preserved.

## Usage

``` r
# S3 method for class 'rdClass'
e1 * e2
```

## Arguments

- e1:

  A scalar numeric value or a `rdClass` object.

- e2:

  A `rdClass` object or a scalar numeric value.

## Value

A `rdClass` object scaled by the numeric scalar.

## Examples

``` r
rd <- rdClass(matrix(1:12, nrow = 4))
rd_scaled <- rd * 0.5
```
