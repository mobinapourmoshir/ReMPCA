# Example for Functional Data (fd)
fd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
fd_object <- fdClass(data = fd_data,
                    argval = NULL,  # Grid points for columns
                    Smoothing_parameter = NULL,  # Custom smoothing parameter
                    Sparsity_parameter = NULL)  # Custom sparsity parameter

# Display the created fd object
print(fd_object)
print(attr(fd_object, "GridPoints_v"))  # Display grid points for columns
print(attr(fd_object, "Smoothing_parameter"))  # Display smoothing parameter

# Example for Regular Data (rd)
rd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # Example regular data matrix (10 rows, 10 columns)
rd_object <- rdClass(data = rd_data,
                    Sparsity_parameter = NULL)  # Custom sparsity parameter

# Display the created rd object
print(rd_object)
print(attr(rd_object, "Sparsity_parameter"))  # Display sparsity parameter

# Example for Hybrid Data (hd)
fd_object2 <- fdClass(data = matrix(rnorm(100), nrow = 10, ncol = 10))  # Another fd object
rd_object2 <- rdClass(data = matrix(rnorm(100), nrow = 10, ncol = 10))  # Another rd object

hd_list <- list(fd_object, rd_object)  # List of fd and rd objects
object_list <- hdClass(hdlist = hd_list,
                     argval = NULL,  # Grid points for rows
                     Smoothing_parameter = NULL,  # Custom smoothing parameter for rows
                     Sparsity_parameter = NULL)  # Custom sparsity parameter for rows

# Display the created hd object
print(object_list)
print(attr(object_list, "GridPoints_u"))  # Display grid points for rows
print(attr(object_list, "Smoothing_parameter"))  # Display row smoothing parameter
print(attr(object_list, "Sparsity_parameter"))  # Display row sparsity parameter

