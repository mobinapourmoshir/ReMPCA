# Indexing operator for imgClass

Enables subsetting of an `imgClass` object by rows.

## Usage

``` r
# S3 method for class 'imgClass'
x[i = NULL]
```

## Arguments

- x:

  An object of class `imgClass`.

- i:

  Row indices (Images). If `NULL`, all images are included.

## Value

A new `imgClass` object with subsetted data and inherited attributes.
