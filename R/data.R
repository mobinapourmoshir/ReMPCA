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






##### Example for running ReMPCA ####
set.seed(123)
# Common left singular vector u (1 x 50)
x <- seq(0,2*pi, len = 150);u <- sin(x); u[1:75] <- 0 #u <- u / sqrt(sum(uˆ2))
# Variable 1
m1 <- 20 ;v1 <- rnorm(m1, sd = 0.8) ;zero_indices_v1 <- sample(1:m1, 5) # Select 5 indices
v1[zero_indices_v1] <- v1[zero_indices_v1] * 1e-2 + rnorm(5, sd = 0.01)
X1 <- outer(u, v1) + rnorm(length(outer(u, v1)), sd = 0.03)

# Variable 2
x2 <- seq(0, 2*pi, length.out = 40) ;v2 <- cos(2*x2) ;v2[20:40] <- 0
X2 <- outer(u, v2) ;X2 <- X2 + rnorm(length(X2), sd = 0.5)

# Variable 3
m3 <- 30; x3 <- seq(0, pi, length.out = m3); v3 <- sin(x3)
X3 <- outer(u, v3) + rnorm(length(outer(u, v3)), sd = 0.2)

# Variable 4
m4 <- 10 ;v4 <- rnorm(m4)
X4 <- outer(u, v4) + rnorm(length(outer(u, v4)), sd = 0.05)


X <- cbind(outer(u, v1), outer(u, v2), outer(u, v3), outer(u, v4))

# Object
rd_object1 <- rdClass(data = as.matrix(X1), Sparsity_parameter = round(seq(1,19, length.out = 10)))
fd_object2 <- fdClass(data = as.matrix(X2), argval = NULL, Smoothing_parameter = NULL, Sparsity_parameter = round(seq(0, 39, length.out = 20)))
fd_object3 <- fdClass(data = as.matrix(X3), argval = NULL, Smoothing_parameter = NULL, Sparsity_parameter = 0) #round(seq(0,29, length.out = 15)))
rd_object4 <- rdClass(data = as.matrix(X4), Sparsity_parameter = 0) #round(seq(0,9, length.out = 5)))
hd_list <- list(rd_object1, fd_object2, fd_object3, rd_object4)
object_list <- hdClass(hdlist = hd_list, argval = NULL, Smoothing_parameter = NULL, Sparsity_parameter = round(seq(0,149, length.out = 20)))

# Arguments
hd = object_list
centerhds = FALSE
num_pcs = 1
smoothness_type = "Second_order"
sparse_tuning_type = "soft"
nfolds_u = 5
nfolds_v = NULL
thresh = 1e-10
maxit = 100
tuning_iter = 1
parallel = FALSE
tuning_order = "Sparsity"
cv.pick = "1se"
sparse_tuning_u = NULL
sparse_tuning_v = NULL
smooth_tuning_u = NULL
smooth_tuning_v = NULL


# Plots
par(mfrow = c(2, 4), pin = c(1, 1), mar = c(2, 2, 2, 2), pty = "s")
matplot(v1, type = 'o', main = 'v1, Regular') ;abline(h = 0, col = 'red', lty = 2)
image(t(X1), main = 'X1, Regular')
matplot(svd(X2)$v[,1], type = 'l', main = 'v2, Functional')
image(t(X2), main = 'X2, Functional')
matplot(v3, type = 'l', main = 'v3, Functional')
image(t(X3), main = 'X3, Functional')
matplot(v4, type = 'o', main = 'v4, Regular')
image(t(X4), main = 'X4, Regular')
