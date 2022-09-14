\name{PXR}
\alias{PXR}
\alias{precision.at.recall.level}
\alias{precision.at.recall.level.over.classes}
\alias{precision.at.multiple.recall.level}
\alias{precision.at.multiple.recall.level.over.classes}
\alias{precision.at.all.recall.levels}

\title{
Precision at a given recall level measures
}
\description{
Set of functions to compute the precision at fixed recall levels.
}
\usage{
precision.at.recall.level(scores, labels, rec.level = 0.2)
precision.at.recall.level.over.classes(target, predicted,  
                          g, rec.level = 0.2, root = "00")
precision.at.multiple.recall.level(scores, labels, 
           rec.levels = seq(from = 0.1, to = 1, by = 0.1))
precision.at.multiple.recall.level.over.classes(target, 
   predicted, rec.levels = seq(from = 0.1, to = 1, by = 0.1))
precision.at.all.recall.levels(scores, labels, resolution=1)
}

\arguments{
  \item{scores}{
 vector of the predicted scores in [0,1]
}
  \item{labels}{
0/1 vector of the true labels 
}
  \item{rec.level}{
rec.level: the desired recall level (def: 0.2)
}
  \item{target}{
matrix with the target multilabels; rows correspond to examples, columns to classes
}
  \item{predicted}{
matrix with the predicted multilabels; rows correspond to examples, columns to classes
}
  \item{g}{
graph of the classes (object of class graphNEL, package graph). 
If missing, no per level results are computed.
}
  \item{root}{
the name of the root node (def. "00") of the graph g.
}
  \item{rec.levels}{
a vector with the desired recall levels (def. 0.1 to 1 by 0.1 step)
}
 \item{resolution}{
a number between 0 and 1 (def. 1). This represents the fraction of precision, recall and f-score values returned. 
}

}
\details{
\code{precision.at.recall.level} computes the precision at a given recall level for a single class.

\code{precision.at.recall.level.over.classes} computes precision at a given recall level for a set of classes.

\code{precision.at.multiple.recall.level} computes the precision at multiple levels of recall for a single class.

\code{precision.at.multiple.recall.level.over.classes} computes the precision at multiple levels of recall for multiple classes.

\code{precision.at.all.recall.levels} compute the precision at all recall levels  for a single class.
It returns a pair of precision and recall values by moving a threshold from the lowest to the highest score: a number of precision and recall values equal to the number n of available examples is returned
if resolution=1, otherwise a number of values equal to n * resolution.
}
\value{
\code{precision.at.recall.level} returns the precision at the requested recall

\code{precision.at.recall.level.over.classes} a list with three elements:
\item{- average }{ the average precision at a given recall level across classes. }
\item{- per.level }{ a named vector with average  precision at a given recall level for each level of the hierarchy; names correspond to levels }
\item{- per.class }{ a named vector with precision at a given recall level for each class. Names correspond to classes}

\code{precision.at.multiple.recall.level} a list with 2 elements:
\item{- precisions }{a vector with the precision at different recall levels}
\item{- f.score }{a vector with the f-score at different recall levels}

\code{precision.at.multiple.recall.level.over.classes} 
\item{- PXR }{a matrix with the precisions at different recall levels: rows are classes, columns    precisions at different recall levels}
\item{- avgPXR }{a vector with the the average precisions at different recall levels across classes}

\code{precision.at.all.recall.levels} a list with 3 elements:
\item{- precision }{precision at different thresholds}
\item{- recall }{recall at different thresholds}
\item{- f.score }{f.score at different thresholds}

}

\seealso{
\code{\link{AUC.measures}}, \code{\link{F.measures}}
}
\examples{
# preparing pseudo-random predictions and target-labels for examples: 
# 100 examples and 10 classes
Scores <- matrix(runif(1000),nrow=100);
Targets <- matrix(integer(1000),nrow=100);
Targets[Scores>0.5] <- 1;
# adding noise to scores
Scores <- Scores + matrix(rnorm(1000, sd=0.3),nrow=100);
colnames(Scores) <-colnames(Targets) <- LETTERS[1:10];
# getting scores and labels of class "A"
scores <- Scores[,"A"];
labels <- Targets[,"A"];
# precsion at 0.4 recall level for class A
precision.at.recall.level(scores, labels, rec.level=0.4);
# precision at 0.4 recall level for all the 10 classes
precision.at.recall.level.over.classes(Targets, Scores, rec.level=0.4);
# precision at multiple recall levels for class A
levels <- seq(from=0.1, to=1, by=0.1);
precision.at.multiple.recall.level(scores, labels, rec.levels=levels);
# precision at multiple recall levels for all the 10 classes
precision.at.multiple.recall.level.over.classes(Targets, Scores);
# precision, recall and f-score for a single class obtained 
# by moving the threshold across the examples
precision.at.all.recall.levels(scores, labels);

}
