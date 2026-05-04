# Indexing Operator for hdClass

Enables subsetting of an `hdClass` object by variables.

## Usage

``` r
# S3 method for class 'hdClass'
x[k = NULL]
```

## Arguments

- x:

  An object of class `hdClass`.

- k:

  Variables indices. If `NULL`, all variables are included.

## Value

A new `hdClass` object with subsetted data and inherited attributes.
