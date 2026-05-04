# Indexing operator for rdClass

Enables subsetting of an `rdClass` object by rows and columns.

## Usage

``` r
# S3 method for class 'rdClass'
x[i = NULL, j = NULL]
```

## Arguments

- x:

  An object of class `rdClass`.

- i:

  Row indices (observations). If `NULL`, all rows are included.

- j:

  Column indices. If `NULL`, all columns are included.

## Value

A new `rdClass` object with subsetted data and inherited attributes.
