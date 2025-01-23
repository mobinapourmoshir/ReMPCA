library(spls)
data(yeast)
F <- yeast$y
S <- yeast$x

matplot(t(F[1:5,]), type = 'l')

F1 <- F[1:5, 1:10]
F2 <- F[1:5, 11:18]
S <- S[1:5,1:3]
hd_obj_data <- hd(fd_matrices = list(F1,F2), nfd_matrices = list(S))
