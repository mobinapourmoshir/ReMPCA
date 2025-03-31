library(spls)
data(yeast)
F <- yeast$y
S <- yeast$x

#matplot(t(F[1:5,]), type = 'l')

F1 <- F[1:5, 1:10]
F2 <- F[1:5, 11:18]
S <- S[1:5,1:3]




# Simulated sample data
### 1. Set basic parameters
d = 100                   # Length of each sin signal
t = seq(0, 1, length.out = d)  # Generate 'd' points from 0 to 1
nn = 300                  # Total number of observations (e.g., rows)
t_total = d + d           # Here t_total = 200, i.e., two segments each of length d

### 2. Construct two "segmented" sin signals v11, v12
#    v11: first d points are sin(pi * t), last d points are sin(pi * 0)=0
#    v12: first d points are 0, last d points are sin(2*pi * t)
v11 = sin(1 * pi * c(t, rep(0, d)))
v12 = sin(2 * pi * c(rep(0, d), t))

### 3. Quick plots to inspect signals
par(mfrow = c(2, 2))
plot(v11, type = "l", main = "v11")
plot(v12, type = "l", main = "v12")


### 4. Construct matrix U1 (with 'nn' rows and 2 columns)
#    We split the 60% of 'nn' rows into three groups, each with different means.
#    U111, U121, U131 -> the first column; U112, U122, U132 -> the second column
U111 <- matrix(rnorm(nn/3, mean = 2, sd = 0.1), ncol=1)
U121 <- matrix(rnorm(nn/3, mean = 1, sd = 0.01), ncol=1)
U131 <- matrix(rnorm(nn/3, mean = 0, sd = 0.1), ncol=1)

U112 <- matrix(rnorm(nn/3, mean = 0, sd = 0.1), ncol=1)
U122 <- matrix(rnorm(nn/3, mean = 1, sd = 0.1), ncol=1)
U132 <- matrix(rnorm(nn/3, mean = 2, sd = 0.1), ncol=1)

# Combine them into a (nn) x 2 matrix
U1_obj = rbind(
  cbind(U111),
  cbind(U121),
  cbind(U131)
)

U2_obj = rbind(
  cbind(U112),
  cbind(U122),
  cbind(U132)
)


### 5. Construct the data matrix X1:

X1 = U1_obj %*% 10 %*% rbind(v11) +
  matrix(
    rnorm(t_total * nn, mean = 0, sd = 6),
    nrow = nn
  )

matplot(t(X1), type = 'l')


X2= U2_obj %*% 3 %*% rbind(v12) +
  matrix(
    rnorm(t_total * nn, mean = 0, sd = 6),
    nrow = nn
  )

matplot(t(X2), type = 'l')

fd_object <- fdClass(data = X1,
                     argval = NULL,  # Grid points for columns
                     Smoothing_parameter = NULL,  # Custom smoothing parameter
                     Sparsity_parameter = round(seq(0,200, length.out = 30)))  # Custom sparsity parameter

rd_object <- rdClass(data = X2,
                     Sparsity_parameter = round(seq(0,200, length.out = 30)))


hd_list <- list(fd_object, rd_object)  # List of fd and rd objects
object_list <- hdClass(hdlist = hd_list,
                       argval = NULL,  # Grid points for rows
                       Smoothing_parameter = NULL,  # Custom smoothing parameter for rows
                       Sparsity_parameter = round(seq(0,300, length.out = 30)))  # Custom sparsity parameter for rows

