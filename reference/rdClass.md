# Regular Data Class

The `rdClass` class denotes an element of regular data.

- Data objects are constructed with observations in rows and variables
  in columns.

- The smoothing parameter is automatically set to zero (no smoothing).

- Users can assign sparsity parameters only!

- No grid points are involved in `rdClass` objects.

## Usage

``` r
rdClass(data, Sparsity_parameter = 0)
```

## Arguments

- data:

  A matrix representing the data, with rows indicating observations and
  columns representing variables.

- Sparsity_parameter:

  - A fixed number representing the level of sparsity for columns, or a
    vector of numerical values that will undergo cross-validation (CV)
    to determine the optimal value.

  - For no sparsity, set it to 0.

  - If `NULL`, the sparsity parameter will be tuned automatically.

## Examples

``` r
# Example for Regular Data (rd)
rd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # 10 rows, 10 columns
rd_object <- rdClass(data = rd_data,
                     Sparsity_parameter = 3)  # Custom sparsity parameter

# Display the created rd object
print(rd_object)
#> Regular Data (rdClass) Object
#> -----------------------------------
#> Dimensions: 10 x 10
#> Sparsity Parameter: 3 
#> -----------------------------------
#> First few rows and columns of the data:
#>              V1          V2         V3         V4          V5
#> [1,]  1.1666841 -1.11449293 -0.4257516 -1.4491557  0.33441525
#> [2,] -0.2718239 -0.84205881  0.6346590 -0.7915040 -0.09349076
#> [3,] -0.3137043 -1.50429421 -0.5785249 -0.5044807  0.30406208
#> [4,] -0.2863002 -0.28402646 -0.1691090  0.4018267 -0.47650753
#> [5,] -0.8034403  0.04286904 -1.9192325  0.9713965 -0.24131172
print(attr(rd_object, "Sparsity_parameter"))  # Display sparsity parameter
#> [1] 3
```
