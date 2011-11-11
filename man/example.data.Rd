\name{example.data}
\alias{example.data}
\alias{Scores}
\alias{T}
\alias{PrecRec}
\docType{data}
\title{ Datasets used in the examples of the package}
\description{
  Collection of datasets used in the examples of the package
}
\usage{
data(Scores)
data(T)
data(PrecRec)
}
\details{
The \code{T} data is a named 1901 X 10 matrix whose rows correspondes to yest genes, while columns correspond to 
10 FunCat (Functional Categories) classes. If T_ij = 1 gene i belong to class j, if T_ij = 0 gene i does not belong to class j.
The \code{Scores} data is a named 1901 X 10 matrix representing scores (likelihood) that a given gene belongs to a given class: higher the value higher the likelihood. 
\code{PrecRec} is  a matrix representing precision at 10 different recall values of 7 methods for gene function prediction.
}

