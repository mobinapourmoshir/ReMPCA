# ReMPCA Smooth and Sparse Multivariate Functional Principal Component Analysis

ReMPCA Smooth and Sparse Multivariate Functional Principal Component
Analysis

## Usage

``` r
ReMPCA(
  hd,
  centerhds = TRUE,
  num_pcs = 1,
  nfolds_u = 5,
  nfolds_v = NULL,
  thresh = 1e-10,
  maxit = 100,
  tuning_iter = 1,
  parallel = FALSE,
  weights = NULL,
  smoothness_type = "Second_order",
  sparse_tuning_type = "soft",
  tuning_order = "Sparsity",
  cv.pick = "1se",
  sparse_tuning_u = NULL,
  sparse_tuning_v = NULL,
  smooth_tuning_u = NULL,
  smooth_tuning_v = NULL
)
```

## Arguments

- hd:

  An hdClass object.

- centerhds:

  A logical; if True, it demeans the data before calculating the
  principal components.

- num_pcs:

  An integer. The number of principal components.

- nfolds_u:

  An integer. It's used in cross validation approach for tuning the
  level of sparsity of rows.

- nfolds_v:

  An integer vector of length `p` (the number of variables), where the
  `i`-th element specifies the number of cross-validation folds to use
  for the `i`-th variable. If not provided (`NULL`), a value of 5 will
  be assigned to all variables by default.

- thresh:

  The convergence threshold in power algorithm.

- maxit:

  Maximum number of iterations in power algorithm.

- tuning_iter:

  Integer specifying the number of iterations to perform during the
  tuning process for conditional smoothing and sparsity parameters.

- parallel:

  Logical; if `TRUE`, parallel computation is used to fit models across
  different sparsity parameter values. Users must register a parallel
  backend beforehand using packages such as doParallel, doMC, or
  similar.

- weights:

  Optional numeric vector of scaling weights. If `NULL`, the function
  automatically computes weights based on the inverse square root of the
  average variance of each variable. If set to `0`, no scaling is
  applied. If a numeric vector is provided, its length must match the
  number of variables, and each element is used to scale the
  corresponding variable. This is used to adjust for scale differences
  across variables in the hybrid data object.

- smoothness_type:

  A character string specifying the method used in smoothing u and/or v,
  must be one of "Second_order" (default), "First_order" or "Indicator".

- sparse_tuning_type:

  A character string specifying the sparse calculation method. Must be
  one of "soft" (default), "hard", or "SCAD".

- tuning_order:

  A character string representing the tuning order. If set to
  'Sparsity', sparsity parameters are tuned first, followed by
  smoothness. If set to 'Smoothness', the order is reversed.

- cv.pick:

  A character string specifying the rule used to select the optimal
  tuning parameter during cross-validation for sparsity. If set to
  `'min'`, the tuning parameter corresponding to the minimum
  cross-validation error is chosen. If set to `'1se'` (the default), the
  1-standard-error rule is applied, selecting the most regularized model
  whose error is within one standard error of the minimum.

- sparse_tuning_u:

  Optional. Specifies the sparsity level(s) for rows:

  - A single non-negative integer for fixed sparsity.

  - A numeric vector of candidate values, to be selected via
    cross-validation (CV).

  - Set to `0` for no sparsity.

  - If `NULL`, the function defaults to the `Sparsity_parameter`
    attribute from the `hdClass` object.

- sparse_tuning_v:

  Optional. A list of length `p` (number of variables), where the `i`-th
  element is either a numeric value or a vector specifying candidate
  sparsity levels for the `i`-th variable. If set to `NULL`, the
  function will default to using the `Sparsity_parameter_col` attribute
  from the input object of class `hdClass`.

- smooth_tuning_u:

  Optional. Specifies smoothing parameter for rows:

  - A single number for fixed smoothness.

  - A numeric vector of candidate values, to be selected via generalized
    cross-validation (GCV).

  - Set to `0` for no smoothness

  - If `NULL`, the function defaults to the `Smoothing_parameter`
    attribute from the `hdClass` object.

- smooth_tuning_v:

  Optional. A list of length `p` (number of variables), where the `i`-th
  element is either a numeric value or a vector specifying candidate
  smoothing parameters for the `i`-th variable. If set to `NULL`, the
  function will default to using the `Smoothing_parameter_col` attribute
  from the input object of class `hdClass`.

## Value

ReconstructedData, PCFunctions, PCScores, OptimalAlphaV, OptimalAlphaU,
OptimalGammaV, OptimalGammaU, GCVResultsV, GCVResultsU, CVResultsV,
CVResultsU, VarianceExplained, variable_types
