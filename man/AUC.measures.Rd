\name{AUC.measures}
\alias{AUC.measures}
\alias{AUC.single}
\alias{AUC.single.over.classes}
\alias{compute.mean.AUC.single.over.classes}
\title{
AUC measures
}
\description{
Set of functions to compute the Area Under the ROC Curve (AUC)
}
\usage{
AUC.single(pred, labels)
AUC.single.over.classes(target, predicted, g, root = "00")
compute.mean.AUC.single.over.classes(y)
}

\arguments{
	\item{pred}{
  numeric vector (scores) of  the values  of the predicted labels
  }
	\item{labels}{
  numeric vector of the true labels (0 negative, 1 positive examples)}
	\item{target}{
  matrix with the target multilabels: rows correspond to examples and columns to classes.
  target[i,j] = 1 if example i belongs to class j, target[i,j] = 0 otherwise.
  }
	\item{predicted}{
  a numeric matrix with predicted values (scores): rows correspond to examples and columns to classes.
  }
	\item{g}{
  a graph of class \emph{graphNEL} (package \pkg{graph}) of the classes. If g is missing no per.level results are computed
  }
	\item{root}{
  the name of the root node (def. "00")
  }
   \item{y}{
    a list of lists. The components of the outer list is a list returned from the function \code{AUC.single.over.classes}
  }
}

\details{
  \code{AUC.single} computes the AUC for a single class.
  
  \code{AUC.single.over.classes} computes AUC for a set of classes, including their average values across classes and
  the average values across the levels of the hierarchy (if any); level 1 classes are at distance 1 from the root,
  level 2 the second level, till to last level correponding to the leaves. Note that if the argument g is missing no   per-level values are computed.
  
  \code{compute.mean.AUC.single.over.classes} compute means across folds of AUC.single.over.classes. It can be used to automatically computed average values (for each class, level, or average across classes) across folds.
}
\value{
  \code{AUC.single} returns a numeric value corresponding to the AUC.
  
  \code{AUC.single.over.classes} returns  a list with three elements:
	\item{- average }{the average AUC across classes}
	\item{- per.level}{a named vector with average  AUC for each level of the hierarchy; names correspond to levels}
	\item{- per.class }{a named vector with AUC for each class; names correspond to classes}
	
  \code{compute.mean.AUC.single.over.classes} returns a list obtained by averaging the results across folds of the input y. 
     The components are:
	\item{- average }{the average AUC across classes}
	\item{- per.level}{a named vector with average  AUC for each level of the hierarchy; names correspond to levels}
	\item{- per.class }{a named vector with AUC for each class; names correspond to classes}
}



\seealso{
\code{\link{F.measures}}, \code{\link{PXR}}
}
\examples{
# preparing pseudo.random scores and target-labels for examples: 100 examples
# and 10 classes
Scores <- matrix(runif(1000),nrow=100);
Targets <- matrix(integer(1000),nrow=100);
Targets[Scores>0.5] <- 1;
# adding noise to scores
Scores <- Scores + matrix(rnorm(1000, sd=0.3),nrow=100);
colnames(Scores) <-colnames(Targets) <- LETTERS[1:10];
# getting scores and labels of class "A"
scores <- Scores[,"A"];
labels <- Targets[,"A"];
# AUC for a single class
AUC.single(scores,labels);
# AUC for the 10 classes
AUC.single.over.classes(Targets, Scores);
}
