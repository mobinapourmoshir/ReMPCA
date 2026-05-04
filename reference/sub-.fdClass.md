# Indexing operator for fdClass

Enables subsetting of an `fdClass` object by rows and columns.

## Usage

``` r
# S3 method for class 'fdClass'
x[i = NULL, j = NULL]
```

## Arguments

- x:

  An object of class `fdClass`.

- i:

  Row indices (observations). If `NULL`, all rows are included.

- j:

  Column indices (grid points). If `NULL`, all columns are included.

## Value

A new `fdClass` object with subsetted data and inherited attributes.
