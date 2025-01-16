library(spls)
data(yeast)
F <- yeast$y
S <- yeast$x


matplot(t(F[1:5,]), type = 'l')


hd_obj_data <- hd(fd_matrices = list(F), nfd_matrices = list(S))

