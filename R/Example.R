# Example for Functional Data (fd)
fd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
fd_object <- fdClaa(fdmatrix = fd_data,
                    argval = NULL,  # Grid points for columns
                    Smoothing_parameter = NULL,  # Custom smoothing parameter
                    Sparsity_parameter = NULL)  # Custom sparsity parameter

# Display the created fd object
print(fd_object)
print(attr(fd_object, "GridPoints_v"))  # Display grid points for columns
print(attr(fd_object, "Smoothing_parameter"))  # Display smoothing parameter

# Example for Regular Data (rd)
rd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # Example regular data matrix (10 rows, 10 columns)
rd_object <- rdClaa(data = rd_data,
                    Sparsity_parameter = 3)  # Custom sparsity parameter

# Display the created rd object
print(rd_object)
print(attr(rd_object, "Sparsity_parameter"))  # Display sparsity parameter

# Example for Hybrid Data (hd)
fd_object2 <- fdClaa(data = matrix(rnorm(100), nrow = 10, ncol = 10))  # Another fd object
rd_object2 <- rdClaa(data = matrix(rnorm(100), nrow = 10, ncol = 10))  # Another rd object

hd_list <- list(fd_object, rd_object)  # List of fd and rd objects
hd_object <- hdClass(hdlist = hd_list,
                     row_argval = seq(0, 1, length.out = 10),  # Grid points for rows
                     row_smoothing_parameter = 0.5,  # Custom smoothing parameter for rows
                     row_sparsity_parameter = 2)  # Custom sparsity parameter for rows

# Display the created hd object
print(hd_object)
print(attr(hd_object, "GridPoints_u"))  # Display grid points for rows
print(attr(hd_object, "row_smoothing_parameter"))  # Display row smoothing parameter
print(attr(hd_object, "row_sparsity_parameter"))  # Display row sparsity parameter
