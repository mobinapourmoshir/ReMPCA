# Compute Scaling Weights for `hdClass` Object

This function computes scaling weights for each component (column) in a
`hdClass` object. The goal is to normalize the contributions of
different functional variables by accounting for their scale, using the
inverse of total variance (integrated over domain).

## Usage

``` r
scale_hd(hd_obj)
```

## Arguments

- hd_obj:

  A hybrid data object of class `hdClass`.

## Value

A numeric vector of weights (length equals the number of variables in
the hybrid data).
