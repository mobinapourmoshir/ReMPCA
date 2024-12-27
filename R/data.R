library(spls)
data(yeast)
F <- yeast$y
S <- yeast$x


matplot(S, type = 'l')
matplot(F, type = 'l')
matplot(t(F), type = 'l')
matplot(t(S), type = 'l')


hd_obj_data <- hd(fd_matrices = list(F), nfd_matrices = list(S))

