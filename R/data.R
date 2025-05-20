library(spls)
data(yeast)
F <- yeast$y
S <- yeast$x

#matplot(t(F[1:5,]), type = 'l')

F1 <- F[1:5, 1:10]
F2 <- F[1:5, 11:18]
S <- S[1:5,1:3]


############## Example ##############
source("fdClass.R")
source("rdClass.R")
source("hdClass.R")
source("imgClass.R")
source("fdMethods.R")
source("rdMethods.R")
source("hdMethods.R")

# Example for Functional Data (fd)
fd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # 10 rows, 10 columns
fd_object <- fdClass(data = fd_data,
                     argval = seq(0, 1, length.out = 10),  # Grid points for columns
                     Smoothing_parameter = 0.5,  # Custom smoothing parameter
                     Sparsity_parameter = 2)  # Custom sparsity parameter

# Display the created fd object
print(fd_object)
print(attr(fd_object, "GridPoints_v"))  # Display grid points for columns
print(attr(fd_object, "Smoothing_parameter"))  # Display smoothing parameter

is.fd(fd_object)
is.rd(fd_object)




# Example for Regular Data (rd)
rd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # 10 rows, 10 columns
rd_object <- rdClass(data = rd_data,
                     Sparsity_parameter = 3)  # Custom sparsity parameter

# Display the created rd object
print(rd_object)
print(attr(rd_object, "Sparsity_parameter"))  # Display sparsity parameter

is.rd(rd_object)
is.fd(rd_object)

convert2fd <- as.fdClass(rd_object)
convert2rd <- as.rdClass(fd_object)

convert2fd <- as.fdClass(rd_object, Smoothing_parameter = c(0.1,0.2,0.3))
convert2rd <- as.rdClass(fd_object, Sparsity_parameter = seq(1:9))

# Example for Hybrid Data (hd)
fd_object2 <- fdClass(data = matrix(rnorm(100), nrow = 10, ncol = 10))  # Another fd object
rd_object2 <- rdClass(data = matrix(rnorm(100), nrow = 10, ncol = 10))  # Another rd object

hd_list <- list(fd_object, rd_object)  # List of fd and rd objects
hd_object <- hdClass(hdlist = hd_list,
                     argval = seq(0, 1, length.out = 10),  # Grid points for rows
                     Smoothing_parameter = 0.5,  # Custom smoothing parameter for rows
                     Sparsity_parameter = 2)  # Custom sparsity parameter for rows

# Display the created hd object
print(hd_object)
print(attr(hd_object, "GridPoints_u"))  # Display grid points for rows
print(attr(hd_object, "Smoothing_parameter"))  # Display row smoothing parameter
print(attr(hd_object, "Sparsity_parameter"))  # Display row sparsity parameter


is.hdClass(hd_object)
is.hdClass(fd_object2)
is.hdClass(rd_object2)
is.rd(hd_object)
is.fd(hd_object)
is.rd(rd_object2)

convert2hd <- as.hdClass(hd_object, Smoothing_parameter = c(0.1,0.2,0.3))
convert2hd <- as.hdClass(fd_object, Smoothing_parameter = c(0.1,0.3))
convert2hd <- as.hdClass(rd_object, Smoothing_parameter = seq(0.01,0.1), Sparsity_parameter = seq(0:9))



# Example for images (imgClass)
set.seed(123)
img1 <- matrix(rnorm(64), 8, 8)
img2 <- matrix(rnorm(64), 8, 8)

# one image
fd_img <- imgClass(img1, Smoothing_parameter = 0.5, Sparsity_parameter = 2)

# Multiple images
rd_img <- imgClass(list(img1, img2), Smoothing_parameter = 0, Sparsity_parameter = NULL)

