\name{sracppfull}
\alias{sracppfull}
\title{Compute the sequential rank agreement between k ranked lists}
\usage{
sracppfull(rankMat)
}
\arguments{
  \item{rankMat}{A matrix with k columns corresponding to
  the k ranked lists. Elements of each column are integers
  between 1 and the length of the lists}
}
\value{
A vector of the same length as the rows in rankMat
containing the sequential rank agreement between the lists
for each depth
}
\description{
Computes the sequential rank agreement (number of items
present in all k lists divided by the current rank) for
each rank in the k lists
}
\author{
Claus Ekstrøm <ekstrom@sund.ku.dk>
}

