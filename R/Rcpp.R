library(Rcpp)


sourceCpp(file="sra.cpp")


nitems <- 20
m <- cbind(1:nitems, c(1,3,4,2,5,6, sample(7:nitems)))

sra(m, 6, B=1, cens=c(6,6))
