library(spls)
data(yeast)
F <- yeast$y
S <- yeast$x


matplot(S, type = 'l')
matplot(F, type = 'l')
matplot(t(F), type = 'l')
matplot(t(S), type = 'l')
